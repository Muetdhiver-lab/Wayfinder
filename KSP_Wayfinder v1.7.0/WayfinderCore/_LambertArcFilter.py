"""Coarse real-ephemeris Lambert filter for Tisserand sequence candidates."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import math

import pykep as pk
from pykep.core import DAY2SEC, lambert_problem

from _SequenceScout import SequenceCandidate
from _Trajectory import fast_ejection_from_vinf


def _linspace(lower, upper, count):
    count = int(count)
    if count < 1:
        raise ValueError("sample count must be positive")
    if count == 1:
        return [0.5 * (float(lower) + float(upper))]
    step = (float(upper) - float(lower)) / (count - 1)
    return [float(lower) + index * step for index in range(count)]


def _norm(vector):
    return math.sqrt(sum(float(value) ** 2 for value in vector))


def _subtract(left, right):
    return [float(a) - float(b) for a, b in zip(left, right)]


@dataclass(frozen=True)
class LambertArc1Config:
    """Sampling and dimensionless acceptance policy for the first arc."""

    t0_samples: int = 16
    tof_samples: int = 16
    transfer_tof_factors: tuple[float, float] = (0.55, 2.25)
    resonant_period_factors: tuple[float, float] = (0.8, 2.2)
    direct_tof_factors: tuple[float, float] = (0.55, 2.0)
    maximum_departure_ratio: float = 1.00
    maximum_tisserand_vinf_relative_error: float = 0.20
    vinf_match_weight: float = 0.50
    parking_altitude_m: float = 100000.0

    def __post_init__(self):
        if int(self.t0_samples) < 2:
            raise ValueError("t0_samples must be at least two")
        if int(self.tof_samples) < 2:
            raise ValueError("tof_samples must be at least two")
        for name in (
            "transfer_tof_factors", "resonant_period_factors",
            "direct_tof_factors",
        ):
            lower, upper = (float(value) for value in getattr(self, name))
            if lower <= 0.0 or upper <= lower:
                raise ValueError(name + " must contain increasing positive values")
        if float(self.maximum_departure_ratio) <= 0.0:
            raise ValueError("maximum_departure_ratio must be positive")
        if float(self.maximum_tisserand_vinf_relative_error) < 0.0:
            raise ValueError(
                "maximum_tisserand_vinf_relative_error must be non-negative"
            )
        if float(self.vinf_match_weight) < 0.0:
            raise ValueError("vinf_match_weight must be non-negative")
        if float(self.parking_altitude_m) < 0.0:
            raise ValueError("parking_altitude_m must be non-negative")

    def to_dict(self):
        values = asdict(self)
        for name in (
            "transfer_tof_factors", "resonant_period_factors",
            "direct_tof_factors",
        ):
            values[name] = list(values[name])
        return values


@dataclass(frozen=True)
class LambertArcSolution:
    encounter_body: str
    t0_days: float
    tof_days: float
    departure_vinf_mps: float
    arrival_vinf_mps: float
    ejection_dv_mps: float
    ejection_strategy: str

    def to_dict(self):
        return asdict(self)


@dataclass(frozen=True)
class LambertArc1Assessment:
    candidate: SequenceCandidate
    accepted: bool
    reason: str
    direct_reference_vinf_mps: float
    direct_reference_ejection_dv_mps: float
    solution: LambertArcSolution | None
    departure_ratio: float | None
    tisserand_vinf_relative_error: float | None
    score: float | None
    model: str = "pykep_zero_rev_lambert_arc1_grid"

    def to_dict(self):
        return {
            "candidate": self.candidate.to_dict(),
            "accepted": bool(self.accepted),
            "reason": self.reason,
            "direct_reference_vinf_mps": self.direct_reference_vinf_mps,
            "direct_reference_ejection_dv_mps": (
                self.direct_reference_ejection_dv_mps
            ),
            "solution": self.solution.to_dict() if self.solution else None,
            "departure_ratio": self.departure_ratio,
            "tisserand_vinf_relative_error": (
                self.tisserand_vinf_relative_error
            ),
            "score": self.score,
            "model": self.model,
        }


@dataclass(frozen=True)
class T0BinCandidate:
    """Best feasible first arc for one Tisserand sequence in one T0 bin."""

    candidate: SequenceCandidate
    bin_start_days: float
    bin_end_days: float
    solution: LambertArcSolution
    direct_reference: LambertArcSolution
    ejection_ratio: float
    tisserand_vinf_relative_error: float
    score: float
    model: str = "tisserand_lambert_arc1_t0_bin"

    def to_dict(self):
        return {
            "candidate": self.candidate.to_dict(),
            "bin_start_days": self.bin_start_days,
            "bin_end_days": self.bin_end_days,
            "solution": self.solution.to_dict(),
            "direct_reference": self.direct_reference.to_dict(),
            "ejection_ratio": self.ejection_ratio,
            "tisserand_vinf_relative_error": (
                self.tisserand_vinf_relative_error
            ),
            "score": self.score,
            "model": self.model,
        }


@dataclass(frozen=True)
class T0BinScanResult:
    t0_start_days: float
    t0_end_days: float
    bin_width_days: float
    direct_reference: LambertArcSolution
    candidates: tuple[T0BinCandidate, ...]
    maximum_ejection_ratio: float
    model: str = "tisserand_lambert_arc1_t0_bin_scan"

    def ranked_unique(self, limit=None):
        """Return each sequence once, using its best bin occurrence."""
        best = {}
        for row in self.candidates:
            sequence = row.candidate.sequence
            current = best.get(sequence)
            if current is None or row.score < current.score:
                best[sequence] = row
        rows = sorted(
            best.values(),
            key=lambda row: (
                row.score,
                row.candidate.proxy_cost_mps,
                row.candidate.sequence,
            ),
        )
        if limit is not None:
            rows = rows[:max(0, int(limit))]
        return rows

    def bin_count(self, sequence):
        sequence = tuple(sequence)
        return sum(
            row.candidate.sequence == sequence for row in self.candidates
        )

    def to_dict(self):
        return {
            "t0_start_days": self.t0_start_days,
            "t0_end_days": self.t0_end_days,
            "bin_width_days": self.bin_width_days,
            "direct_reference": self.direct_reference.to_dict(),
            "candidates": [row.to_dict() for row in self.candidates],
            "maximum_ejection_ratio": self.maximum_ejection_ratio,
            "model": self.model,
        }


class LambertArc1Filter:
    """Reject Tisserand sequences with no plausible real first arc."""

    MODEL = "pykep_zero_rev_lambert_arc1_grid"

    @staticmethod
    def _scored_solution_key(item):
        score, ratio, error, solution = item
        return score, ratio, error, solution.t0_days, solution.tof_days

    def __init__(self, planet_pack_module, config=None):
        self.planet_pack = planet_pack_module
        self.config = config or LambertArc1Config()
        self.bodies = dict(planet_pack_module.BODIES)
        self.edy_to_local_days = float(planet_pack_module.EDY_TO_KDY)
        if not self.bodies:
            raise ValueError("planet pack contains no Lambert bodies")
        first = next(iter(self.bodies.values()))
        self.central_mu = float(first.mu_central_body)

    def _body(self, name):
        try:
            return self.bodies[str(name)]
        except KeyError as exc:
            raise ValueError("Unknown Lambert body: " + str(name)) from exc

    def _orbit_radius(self, name):
        return float(self._body(name).elements(0)[0])

    def _hohmann_time_days(self, departure, arrival):
        departure_radius = self._orbit_radius(departure)
        arrival_radius = self._orbit_radius(arrival)
        transfer_axis = 0.5 * (departure_radius + arrival_radius)
        seconds = math.pi * math.sqrt(
            transfer_axis ** 3 / self.central_mu
        )
        return seconds / DAY2SEC * self.edy_to_local_days

    def _period_days(self, body_name):
        radius = self._orbit_radius(body_name)
        seconds = 2.0 * math.pi * math.sqrt(radius ** 3 / self.central_mu)
        return seconds / DAY2SEC * self.edy_to_local_days

    def _tof_bounds(self, departure, arrival, direct_reference=False):
        if departure == arrival:
            reference = self._period_days(departure)
            factors = self.config.resonant_period_factors
        else:
            reference = self._hohmann_time_days(departure, arrival)
            factors = (
                self.config.direct_tof_factors
                if direct_reference
                else self.config.transfer_tof_factors
            )
        return reference * float(factors[0]), reference * float(factors[1])

    def _solve_grid(
        self,
        departure,
        arrival,
        t0_bounds_days,
        tof_bounds_days,
        t0_samples=None,
        tof_samples=None,
    ):
        departure_body = self._body(departure)
        arrival_body = self._body(arrival)
        t0_lower, t0_upper = (float(value) for value in t0_bounds_days)
        tof_lower, tof_upper = (float(value) for value in tof_bounds_days)
        if t0_upper <= t0_lower:
            raise ValueError("t0_bounds_days must be increasing")
        if tof_upper <= tof_lower or tof_lower <= 0.0:
            raise ValueError("tof_bounds_days must be increasing and positive")

        solutions = []
        t0_samples = self.config.t0_samples if t0_samples is None else int(t0_samples)
        tof_samples = self.config.tof_samples if tof_samples is None else int(tof_samples)
        for t0_days in _linspace(t0_lower, t0_upper, t0_samples):
            departure_epoch_days = t0_days / self.edy_to_local_days
            departure_state = departure_body.eph(pk.epoch(departure_epoch_days))
            for tof_days in _linspace(
                tof_lower, tof_upper, tof_samples,
            ):
                tof_earth_days = tof_days / self.edy_to_local_days
                arrival_epoch_days = departure_epoch_days + tof_earth_days
                arrival_state = arrival_body.eph(pk.epoch(arrival_epoch_days))
                try:
                    lambert = lambert_problem(
                        departure_state[0],
                        arrival_state[0],
                        tof_earth_days * DAY2SEC,
                        self.central_mu,
                        cw=False,
                        multi_revs=0,
                    )
                    departure_velocity = lambert.v0[0]
                    arrival_velocity = lambert.v1[0]
                except (RuntimeError, ValueError, OverflowError):
                    continue
                departure_vinf_vector = _subtract(
                    departure_velocity, departure_state[1],
                )
                arrival_vinf_vector = _subtract(
                    arrival_velocity, arrival_state[1],
                )
                ejection = fast_ejection_from_vinf(
                    departure_body,
                    departure_vinf_vector,
                    self.config.parking_altitude_m,
                )
                solutions.append(LambertArcSolution(
                    encounter_body=arrival,
                    t0_days=t0_days,
                    tof_days=tof_days,
                    departure_vinf_mps=_norm(departure_vinf_vector),
                    arrival_vinf_mps=_norm(arrival_vinf_vector),
                    ejection_dv_mps=float(ejection["dv"]),
                    ejection_strategy=str(ejection["strategy"]),
                ))
        return solutions

    def _direct_reference(self, start, target, t0_bounds_days):
        solutions = self._solve_grid(
            start,
            target,
            t0_bounds_days,
            self._tof_bounds(start, target, direct_reference=True),
        )
        if not solutions:
            raise RuntimeError("No direct Lambert reference solution found")
        return min(solutions, key=lambda solution: solution.ejection_dv_mps)

    @staticmethod
    def _first_flyby_vinf(candidate):
        if not candidate.flybys:
            return None
        return float(candidate.flybys[0].incoming_vinf_mps)

    def assess(self, candidates, t0_bounds_days):
        candidates = list(candidates)
        if not candidates:
            return []
        start_names = {candidate.sequence[0] for candidate in candidates}
        target_names = {candidate.sequence[-1] for candidate in candidates}
        if len(start_names) != 1 or len(target_names) != 1:
            raise ValueError(
                "all candidates must share one departure and one target"
            )
        start = next(iter(start_names))
        target = next(iter(target_names))
        direct_reference = self._direct_reference(
            start, target, t0_bounds_days,
        )
        reference_vinf = direct_reference.departure_vinf_mps
        reference_ejection = direct_reference.ejection_dv_mps

        solutions_by_first_body = {}
        for first_body in sorted({
            candidate.sequence[1] for candidate in candidates
        }):
            solutions_by_first_body[first_body] = self._solve_grid(
                start,
                first_body,
                t0_bounds_days,
                self._tof_bounds(start, first_body),
            )

        assessments = []
        for candidate in candidates:
            first_body = candidate.sequence[1]
            desired_vinf = self._first_flyby_vinf(candidate)
            scored = []
            for solution in solutions_by_first_body[first_body]:
                departure_ratio = (
                    solution.ejection_dv_mps / max(reference_ejection, 1e-12)
                )
                if desired_vinf is None:
                    relative_error = 0.0
                else:
                    relative_error = abs(
                        solution.arrival_vinf_mps - desired_vinf
                    ) / max(desired_vinf, 100.0)
                score = (
                    departure_ratio
                    + float(self.config.vinf_match_weight) * relative_error
                )
                scored.append((score, departure_ratio, relative_error, solution))
            if not scored:
                assessments.append(LambertArc1Assessment(
                    candidate=candidate,
                    accepted=False,
                    reason="no_lambert_solution",
                    direct_reference_vinf_mps=reference_vinf,
                    direct_reference_ejection_dv_mps=reference_ejection,
                    solution=None,
                    departure_ratio=None,
                    tisserand_vinf_relative_error=None,
                    score=None,
                ))
                continue
            feasible = [
                item for item in scored
                if item[1] <= self.config.maximum_departure_ratio
                and (
                    desired_vinf is None
                    or item[2]
                    <= self.config.maximum_tisserand_vinf_relative_error
                )
            ]
            if feasible:
                score, departure_ratio, relative_error, solution = min(
                    feasible, key=self._scored_solution_key,
                )
                accepted = True
                reason = "accepted"
            else:
                score, departure_ratio, relative_error, solution = min(
                    scored, key=self._scored_solution_key,
                )
                accepted = False
                if departure_ratio > self.config.maximum_departure_ratio:
                    reason = "departure_energy_above_direct_reference"
                else:
                    reason = "first_flyby_vinf_mismatch"
            assessments.append(LambertArc1Assessment(
                candidate=candidate,
                accepted=accepted,
                reason=reason,
                direct_reference_vinf_mps=reference_vinf,
                direct_reference_ejection_dv_mps=reference_ejection,
                solution=solution,
                departure_ratio=departure_ratio,
                tisserand_vinf_relative_error=relative_error,
                score=score,
            ))
        return sorted(
            assessments,
            key=lambda assessment: (
                not assessment.accepted,
                math.inf if assessment.score is None else assessment.score,
                assessment.candidate.proxy_cost_mps,
                assessment.candidate.sequence,
            ),
        )

    def filter(self, candidates, t0_bounds_days):
        return [
            assessment
            for assessment in self.assess(candidates, t0_bounds_days)
            if assessment.accepted
        ]

    def scan_t0_bins(
        self,
        candidates,
        t0_bounds_days,
        bin_width_days=100.0,
        t0_step_days=10.0,
        tof_samples=40,
        maximum_ejection_ratio=1.05,
        turn_fraction_weight=0.10,
        flyby_count_weight=0.02,
        include_direct=False,
    ):
        """Scan first-arc feasibility independently in fixed T0 bins.

        The direct reference is global over ``t0_bounds_days``. Each candidate
        is retained in a bin only when one Lambert arc matches its first
        Tisserand encounter and its exact parking-orbit ejection remains below
        ``maximum_ejection_ratio`` times that global direct ejection.
        """
        candidates = list(candidates)
        if not candidates:
            raise ValueError("at least one Tisserand candidate is required")
        t0_start, t0_end = (float(value) for value in t0_bounds_days)
        bin_width_days = float(bin_width_days)
        t0_step_days = float(t0_step_days)
        maximum_ejection_ratio = float(maximum_ejection_ratio)
        if t0_end <= t0_start:
            raise ValueError("t0_bounds_days must be increasing")
        if bin_width_days <= 0.0 or t0_step_days <= 0.0:
            raise ValueError("bin width and T0 step must be positive")
        if int(tof_samples) < 2:
            raise ValueError("tof_samples must be at least two")
        if maximum_ejection_ratio <= 0.0:
            raise ValueError("maximum_ejection_ratio must be positive")

        start_names = {candidate.sequence[0] for candidate in candidates}
        target_names = {candidate.sequence[-1] for candidate in candidates}
        if len(start_names) != 1 or len(target_names) != 1:
            raise ValueError(
                "all candidates must share one departure and one target"
            )
        start = next(iter(start_names))
        target = next(iter(target_names))

        global_t0_samples = max(
            2, int(math.ceil((t0_end - t0_start) / t0_step_days)) + 1,
        )
        direct_solutions = self._solve_grid(
            start,
            target,
            [t0_start, t0_end],
            self._tof_bounds(start, target, direct_reference=True),
            t0_samples=global_t0_samples,
            tof_samples=tof_samples,
        )
        if not direct_solutions:
            raise RuntimeError("No direct Lambert reference solution found")
        direct_reference = min(
            direct_solutions, key=lambda solution: solution.ejection_dv_mps,
        )

        rows = []
        bin_start = t0_start
        while bin_start < t0_end - 1e-12:
            bin_end = min(t0_end, bin_start + bin_width_days)
            bin_t0_samples = max(
                2, int(math.ceil((bin_end - bin_start) / t0_step_days)) + 1,
            )
            solutions_by_first_body = {}
            first_bodies = {
                candidate.sequence[1] for candidate in candidates
                if include_direct or len(candidate.sequence) > 2
            }
            for first_body in sorted(first_bodies):
                solutions = self._solve_grid(
                    start,
                    first_body,
                    [bin_start, bin_end],
                    self._tof_bounds(start, first_body),
                    t0_samples=bin_t0_samples,
                    tof_samples=tof_samples,
                )
                # Bins are lower-inclusive and upper-exclusive, except the last.
                if bin_end < t0_end - 1e-12:
                    solutions = [
                        solution for solution in solutions
                        if solution.t0_days < bin_end - 1e-9
                    ]
                solutions_by_first_body[first_body] = solutions

            for candidate in candidates:
                if not include_direct and len(candidate.sequence) <= 2:
                    continue
                desired_vinf = self._first_flyby_vinf(candidate)
                feasible = []
                for solution in solutions_by_first_body[candidate.sequence[1]]:
                    ejection_ratio = (
                        solution.ejection_dv_mps
                        / max(direct_reference.ejection_dv_mps, 1e-12)
                    )
                    if desired_vinf is None:
                        relative_error = 0.0
                    else:
                        relative_error = abs(
                            solution.arrival_vinf_mps - desired_vinf
                        ) / max(desired_vinf, 100.0)
                    if ejection_ratio > maximum_ejection_ratio:
                        continue
                    if (
                        desired_vinf is not None
                        and relative_error
                        > self.config.maximum_tisserand_vinf_relative_error
                    ):
                        continue
                    flyby_count = max(0, len(candidate.sequence) - 2)
                    score = (
                        ejection_ratio
                        + float(self.config.vinf_match_weight) * relative_error
                        + float(turn_fraction_weight) * candidate.max_turn_fraction
                        + float(flyby_count_weight) * flyby_count
                    )
                    feasible.append((
                        score, ejection_ratio, relative_error, solution,
                    ))
                if not feasible:
                    continue
                score, ejection_ratio, relative_error, solution = min(
                    feasible, key=self._scored_solution_key,
                )
                rows.append(T0BinCandidate(
                    candidate=candidate,
                    bin_start_days=bin_start,
                    bin_end_days=bin_end,
                    solution=solution,
                    direct_reference=direct_reference,
                    ejection_ratio=ejection_ratio,
                    tisserand_vinf_relative_error=relative_error,
                    score=score,
                ))
            bin_start = bin_end

        return T0BinScanResult(
            t0_start_days=t0_start,
            t0_end_days=t0_end,
            bin_width_days=bin_width_days,
            direct_reference=direct_reference,
            candidates=tuple(sorted(
                rows,
                key=lambda row: (
                    row.bin_start_days,
                    row.score,
                    row.candidate.sequence,
                ),
            )),
            maximum_ejection_ratio=maximum_ejection_ratio,
        )
