"""Unphased Tisserand/tree scout for preliminary MGA sequence selection.

The scout intentionally uses the circular, coplanar patched-conic model. It
answers a narrow question: can a chain of unpowered flybys connect heliocentric
orbit-energy regions? It does *not* prove that the encounters can be phased.
The candidates are meant to feed a Lambert first-arc filter and then Wayfinder's
normal L0/funnel optimizer.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
import math


TWO_PI = 2.0 * math.pi


def _linspace(lower, upper, count, endpoint=True):
    count = int(count)
    if count < 1:
        raise ValueError("sample count must be positive")
    if count == 1:
        return [0.5 * (float(lower) + float(upper))]
    denominator = count - 1 if endpoint else count
    step = (float(upper) - float(lower)) / denominator
    return [float(lower) + index * step for index in range(count)]


def _wrap_angle(angle):
    return (float(angle) + math.pi) % TWO_PI - math.pi


def tisserand_parameter(
    semi_major_axis, eccentricity, planet_orbit_radius, inclination=0.0,
):
    """Return the classical dimensionless Tisserand parameter.

    ``planet_orbit_radius`` is the circular reference orbit of the encounter
    body. The formula is meaningful only within the restricted circular model.
    """
    semi_major_axis = float(semi_major_axis)
    eccentricity = float(eccentricity)
    planet_orbit_radius = float(planet_orbit_radius)
    if semi_major_axis <= 0.0 or planet_orbit_radius <= 0.0:
        raise ValueError("semi-major axes must be positive")
    angular_momentum_term = (
        semi_major_axis
        * max(0.0, 1.0 - eccentricity * eccentricity)
        / planet_orbit_radius
    )
    return (
        planet_orbit_radius / semi_major_axis
        + 2.0 * math.sqrt(angular_momentum_term) * math.cos(float(inclination))
    )


def maximum_flyby_turn_angle(vinf_mps, body_mu, periapsis_radius):
    """Maximum unpowered hyperbolic turn angle, in radians."""
    vinf_mps = float(vinf_mps)
    body_mu = float(body_mu)
    periapsis_radius = float(periapsis_radius)
    if vinf_mps <= 0.0:
        return math.pi
    if body_mu <= 0.0 or periapsis_radius <= 0.0:
        raise ValueError("body_mu and periapsis_radius must be positive")
    eccentricity = 1.0 + periapsis_radius * vinf_mps * vinf_mps / body_mu
    return 2.0 * math.asin(min(1.0, 1.0 / eccentricity))


def required_flyby_periapsis(vinf_mps, turn_angle, body_mu):
    """Periapsis radius required for an unpowered turn, or infinity at zero."""
    vinf_mps = float(vinf_mps)
    turn_angle = abs(float(turn_angle))
    body_mu = float(body_mu)
    if turn_angle <= 1e-15 or vinf_mps <= 0.0:
        return math.inf
    sine = math.sin(min(math.pi, turn_angle) / 2.0)
    return body_mu / (vinf_mps * vinf_mps) * (1.0 / sine - 1.0)


@dataclass(frozen=True)
class TisserandScoutConfig:
    """Serializable sampling and pruning policy for :class:`SequenceScout`."""

    max_bodies: int = 5
    departure_vinf_samples: int = 13
    departure_direction_samples: int = 72
    flyby_turn_samples: int = 9
    states_per_sequence: int = 128
    max_states_per_depth: int = 50000
    max_visits_per_body: int = 3
    crossing_tolerance: float = 0.002
    radial_overshoot_fraction: float = 0.10
    departure_vinf_min_mps: float | None = None
    departure_vinf_max_mps: float | None = None
    max_encounter_vinf_mps: float = 6000.0
    flyby_clearance_m: float = 0.0
    flyby_count_penalty_mps: float = 50.0
    turn_fraction_penalty_mps: float = 100.0

    def __post_init__(self):
        for name in (
            "max_bodies", "departure_vinf_samples",
            "departure_direction_samples", "flyby_turn_samples",
            "states_per_sequence", "max_states_per_depth",
            "max_visits_per_body",
        ):
            if int(getattr(self, name)) < 1:
                raise ValueError(name + " must be positive")
        if self.max_bodies < 2:
            raise ValueError("max_bodies must be at least two")
        if not 0.0 <= float(self.crossing_tolerance) <= 0.25:
            raise ValueError("crossing_tolerance must be in [0, 0.25]")
        if float(self.radial_overshoot_fraction) < 0.0:
            raise ValueError("radial_overshoot_fraction must be non-negative")
        if float(self.max_encounter_vinf_mps) <= 0.0:
            raise ValueError("max_encounter_vinf_mps must be positive")
        if float(self.flyby_clearance_m) < 0.0:
            raise ValueError("flyby_clearance_m must be non-negative")
        lower = self.departure_vinf_min_mps
        upper = self.departure_vinf_max_mps
        if lower is not None and float(lower) < 0.0:
            raise ValueError("departure_vinf_min_mps must be non-negative")
        if upper is not None and float(upper) <= 0.0:
            raise ValueError("departure_vinf_max_mps must be positive")
        if lower is not None and upper is not None and float(lower) >= float(upper):
            raise ValueError("departure v-infinity bounds must be increasing")

    def to_dict(self):
        return asdict(self)


@dataclass(frozen=True)
class FlybyEvidence:
    body: str
    incoming_vinf_mps: float
    incoming_angle_deg: float
    selected_turn_deg: float
    maximum_turn_deg: float
    turn_fraction: float
    required_periapsis_radius_m: float | None
    safe_periapsis_radius_m: float
    operational_periapsis_radius_m: float
    tisserand_parameter: float

    def to_dict(self):
        return asdict(self)


@dataclass(frozen=True)
class SequenceCandidate:
    sequence: tuple[str, ...]
    proxy_cost_mps: float
    departure_vinf_mps: float
    arrival_vinf_mps: float
    max_turn_fraction: float
    terminal_periapsis_m: float
    terminal_apoapsis_m: float
    flybys: tuple[FlybyEvidence, ...] = field(default_factory=tuple)
    model: str = "tisserand_circular_coplanar_unphased"

    def to_dict(self):
        values = asdict(self)
        values["sequence"] = list(self.sequence)
        values["flybys"] = [flyby.to_dict() for flyby in self.flybys]
        return values


@dataclass(frozen=True)
class _BodyData:
    name: str
    orbit_radius: float
    orbital_speed: float
    mu_self: float
    safe_radius: float


@dataclass(frozen=True)
class _OrbitState:
    sequence: tuple[str, ...]
    body: str
    semi_major_axis: float
    eccentricity: float
    periapsis: float
    apoapsis: float
    angular_momentum: float
    radial_velocity: float
    vinf_mps: float
    relative_angle: float
    departure_vinf_mps: float
    max_turn_fraction: float
    flybys: tuple[FlybyEvidence, ...]


class SequenceScout:
    """Search unphased MGA body sequences using a sampled Tisserand tree."""

    MODEL = "tisserand_circular_coplanar_unphased"

    def __init__(self, planet_pack_module, config=None):
        self.planet_pack = planet_pack_module
        self.config = config or TisserandScoutConfig()
        bodies = dict(planet_pack_module.BODIES)
        if not bodies:
            raise ValueError("planet pack contains no scoutable bodies")
        first = next(iter(bodies.values()))
        self.central_mu = float(first.mu_central_body)
        self.bodies = {}
        for name, body in bodies.items():
            radius = float(body.elements(0)[0])
            if not math.isclose(
                float(body.mu_central_body), self.central_mu, rel_tol=1e-12,
            ):
                continue
            self.bodies[name] = _BodyData(
                name=name,
                orbit_radius=radius,
                orbital_speed=math.sqrt(self.central_mu / radius),
                mu_self=float(body.mu_self),
                safe_radius=float(body.safe_radius),
            )

    def _body(self, name):
        try:
            return self.bodies[str(name)]
        except KeyError as exc:
            raise ValueError("Unknown scout body: " + str(name)) from exc

    def _default_departure_vinf_bounds(self, start, target):
        start_body = self._body(start)
        target_body = self._body(target)
        transfer_axis = 0.5 * (
            start_body.orbit_radius + target_body.orbit_radius
        )
        transfer_speed = math.sqrt(
            self.central_mu
            * (2.0 / start_body.orbit_radius - 1.0 / transfer_axis)
        )
        direct_vinf = abs(transfer_speed - start_body.orbital_speed)
        lower = max(100.0, 0.25 * direct_vinf)
        upper = max(lower + 100.0, 1.35 * direct_vinf)
        return lower, upper

    def _departure_vinf_bounds(self, start, target):
        default_lower, default_upper = self._default_departure_vinf_bounds(
            start, target,
        )
        lower = self.config.departure_vinf_min_mps
        upper = self.config.departure_vinf_max_mps
        return (
            default_lower if lower is None else float(lower),
            default_upper if upper is None else float(upper),
        )

    def _candidate_body_names(self, start, target, candidate_bodies):
        start_radius = self._body(start).orbit_radius
        target_radius = self._body(target).orbit_radius
        radial_lower = min(start_radius, target_radius)
        radial_upper = max(start_radius, target_radius)
        overshoot = float(self.config.radial_overshoot_fraction)
        if candidate_bodies is None:
            names = [
                name for name, body in self.bodies.items()
                if radial_lower * (1.0 - overshoot)
                <= body.orbit_radius
                <= radial_upper * (1.0 + overshoot)
            ]
            # A first assist on the opposite side of the departure orbit is a
            # common way to change energy (for example Kerbin-Eve-...-Jool).
            # Include the nearest such body without opening the tree to every
            # radial overshoot in the planet pack.
            if target_radius > start_radius:
                opposite = [
                    body for body in self.bodies.values()
                    if body.orbit_radius < start_radius
                ]
                if opposite:
                    names.append(max(opposite, key=lambda body: body.orbit_radius).name)
            elif target_radius < start_radius:
                opposite = [
                    body for body in self.bodies.values()
                    if body.orbit_radius > start_radius
                ]
                if opposite:
                    names.append(min(opposite, key=lambda body: body.orbit_radius).name)
        else:
            names = [str(name) for name in candidate_bodies]
            for name in names:
                self._body(name)
        for required in (start, target):
            if required not in names:
                names.append(required)
        return sorted(set(names), key=lambda name: self._body(name).orbit_radius)

    def _orbit_from_relative_state(
        self, body_name, vinf_mps, relative_angle, sequence,
        departure_vinf_mps, max_turn_fraction, flybys,
    ):
        body = self._body(body_name)
        vinf_mps = float(vinf_mps)
        relative_angle = _wrap_angle(relative_angle)
        tangential_velocity = (
            body.orbital_speed + vinf_mps * math.cos(relative_angle)
        )
        radial_velocity = vinf_mps * math.sin(relative_angle)
        energy = 0.5 * (
            tangential_velocity * tangential_velocity
            + radial_velocity * radial_velocity
        ) - self.central_mu / body.orbit_radius
        if energy >= 0.0:
            return None
        semi_major_axis = -self.central_mu / (2.0 * energy)
        angular_momentum = body.orbit_radius * tangential_velocity
        eccentricity_squared = (
            1.0
            + 2.0 * energy * angular_momentum * angular_momentum
            / (self.central_mu * self.central_mu)
        )
        if eccentricity_squared < -1e-10:
            return None
        eccentricity = math.sqrt(max(0.0, eccentricity_squared))
        return _OrbitState(
            sequence=tuple(sequence),
            body=body_name,
            semi_major_axis=semi_major_axis,
            eccentricity=eccentricity,
            periapsis=semi_major_axis * (1.0 - eccentricity),
            apoapsis=semi_major_axis * (1.0 + eccentricity),
            angular_momentum=angular_momentum,
            radial_velocity=radial_velocity,
            vinf_mps=vinf_mps,
            relative_angle=relative_angle,
            departure_vinf_mps=float(departure_vinf_mps),
            max_turn_fraction=float(max_turn_fraction),
            flybys=tuple(flybys),
        )

    def _encounter_states(self, state, next_body_name):
        body = self._body(next_body_name)
        tolerance = float(self.config.crossing_tolerance)
        radius = body.orbit_radius
        if (
            radius < state.periapsis * (1.0 - tolerance)
            or radius > state.apoapsis * (1.0 + tolerance)
        ):
            return []
        speed_squared = self.central_mu * (
            2.0 / radius - 1.0 / state.semi_major_axis
        )
        tangential_velocity = state.angular_momentum / radius
        radial_squared = speed_squared - tangential_velocity * tangential_velocity
        numerical_scale = max(speed_squared, 1.0)
        if radial_squared < -1e-10 * numerical_scale:
            return []
        radial_speed = math.sqrt(max(0.0, radial_squared))
        signs = (1.0, -1.0) if radial_speed > 1e-12 else (0.0,)
        encounters = []
        seen = set()
        for sign in signs:
            radial_velocity = sign * radial_speed
            relative_tangential = tangential_velocity - body.orbital_speed
            vinf = math.hypot(relative_tangential, radial_velocity)
            angle = math.atan2(radial_velocity, relative_tangential)
            key = (round(vinf, 9), round(angle, 12))
            if key not in seen:
                encounters.append((vinf, angle))
                seen.add(key)
        return encounters

    def _flyby_evidence(self, state, body_name, vinf, incoming_angle, turn):
        body = self._body(body_name)
        maximum_turn = maximum_flyby_turn_angle(
            vinf,
            body.mu_self,
            body.safe_radius + float(self.config.flyby_clearance_m),
        )
        turn_fraction = (
            abs(turn) / maximum_turn if maximum_turn > 1e-15 else 0.0
        )
        parameter = 3.0 - (vinf / body.orbital_speed) ** 2
        required_periapsis = required_flyby_periapsis(
            vinf, turn, body.mu_self,
        )
        return FlybyEvidence(
            body=body_name,
            incoming_vinf_mps=float(vinf),
            incoming_angle_deg=math.degrees(incoming_angle),
            selected_turn_deg=math.degrees(turn),
            maximum_turn_deg=math.degrees(maximum_turn),
            turn_fraction=float(turn_fraction),
            required_periapsis_radius_m=(
                required_periapsis
                if math.isfinite(required_periapsis)
                else None
            ),
            safe_periapsis_radius_m=body.safe_radius,
            operational_periapsis_radius_m=(
                body.safe_radius + float(self.config.flyby_clearance_m)
            ),
            tisserand_parameter=float(parameter),
        )

    def _state_merit(self, state):
        flyby_count = max(0, len(state.sequence) - 1)
        return (
            state.departure_vinf_mps
            + float(self.config.flyby_count_penalty_mps) * flyby_count
            + float(self.config.turn_fraction_penalty_mps)
            * state.max_turn_fraction
        )

    def _state_grid_key(self, state):
        body_radius = self._body(state.body).orbit_radius
        return (
            round(math.log(state.semi_major_axis / body_radius) / 0.08),
            round(state.eccentricity / 0.04),
            round((_wrap_angle(state.relative_angle) + math.pi) / math.radians(8.0)),
            round(state.vinf_mps / 100.0),
        )

    def _select_diverse(self, states, limit):
        """Keep low-cost representatives plus a deterministic grid spread."""
        representatives = {}
        for state in states:
            key = self._state_grid_key(state)
            merit = self._state_merit(state)
            current = representatives.get(key)
            if current is None or merit < current[0]:
                representatives[key] = (merit, state)
        keyed = sorted(representatives.items())
        if len(keyed) <= limit:
            return [entry[1][1] for entry in keyed]

        best_count = max(1, int(limit) // 4)
        selected = [
            state for _, state in sorted(
                representatives.values(), key=lambda item: item[0],
            )[:best_count]
        ]
        spread_count = int(limit) - len(selected)
        pool = [entry[1][1] for entry in keyed]
        if spread_count > 0:
            indices = _linspace(0, len(pool) - 1, spread_count)
            selected.extend(pool[int(round(index))] for index in indices)

        unique = []
        seen = set()
        for state in selected:
            identity = id(state)
            if identity not in seen:
                unique.append(state)
                seen.add(identity)
        if len(unique) < limit:
            for state in pool:
                identity = id(state)
                if identity not in seen:
                    unique.append(state)
                    seen.add(identity)
                if len(unique) >= limit:
                    break
        return unique[:int(limit)]

    def _prune_frontier(self, states):
        by_sequence = {}
        for state in states:
            by_sequence.setdefault(state.sequence, []).append(state)
        selected_by_sequence = {
            sequence: self._select_diverse(
                sequence_states, self.config.states_per_sequence,
            )
            for sequence, sequence_states in sorted(by_sequence.items())
        }
        total = sum(len(states) for states in selected_by_sequence.values())
        if total <= self.config.max_states_per_depth:
            return [
                state for states in selected_by_sequence.values() for state in states
            ]

        # Round-robin truncation preserves every sequence before adding depth.
        frontier = []
        level = 0
        while len(frontier) < self.config.max_states_per_depth:
            added = False
            for sequence in sorted(selected_by_sequence):
                sequence_states = selected_by_sequence[sequence]
                if level < len(sequence_states):
                    frontier.append(sequence_states[level])
                    added = True
                    if len(frontier) >= self.config.max_states_per_depth:
                        break
            if not added:
                break
            level += 1
        return frontier

    def _candidate_from_encounter(self, state, target, vinf):
        flyby_count = max(0, len(state.sequence) - 1)
        proxy = (
            state.departure_vinf_mps
            + float(self.config.flyby_count_penalty_mps) * flyby_count
            + float(self.config.turn_fraction_penalty_mps)
            * state.max_turn_fraction
        )
        return SequenceCandidate(
            sequence=state.sequence + (target,),
            proxy_cost_mps=float(proxy),
            departure_vinf_mps=state.departure_vinf_mps,
            arrival_vinf_mps=float(vinf),
            max_turn_fraction=state.max_turn_fraction,
            terminal_periapsis_m=state.periapsis,
            terminal_apoapsis_m=state.apoapsis,
            flybys=state.flybys,
        )

    def scout(self, start, target, candidate_bodies=None, limit=None):
        """Return unique energy-plausible sequences sorted by proxy cost.

        The result ignores epoch, time of flight, planetary phase and orbital
        inclination. A candidate must therefore pass the future Lambert arc-1
        filter before it is submitted to L0.
        """
        start = str(start)
        target = str(target)
        if start == target:
            raise ValueError("start and target must differ")
        body_names = self._candidate_body_names(
            start, target, candidate_bodies,
        )
        lower_vinf, upper_vinf = self._departure_vinf_bounds(start, target)
        if lower_vinf >= upper_vinf:
            raise ValueError("departure v-infinity bounds must be increasing")

        frontier = []
        directions = _linspace(
            -math.pi, math.pi,
            self.config.departure_direction_samples,
            endpoint=False,
        )
        for vinf in _linspace(
            lower_vinf, upper_vinf, self.config.departure_vinf_samples,
        ):
            for angle in directions:
                state = self._orbit_from_relative_state(
                    start, vinf, angle, (start,), vinf, 0.0, (),
                )
                if state is not None:
                    frontier.append(state)

        best_by_sequence = {}
        for _next_body_count in range(2, self.config.max_bodies + 1):
            next_frontier = []
            for state in frontier:
                for next_body in body_names:
                    if state.sequence.count(next_body) >= self.config.max_visits_per_body:
                        continue
                    for vinf, incoming_angle in self._encounter_states(
                        state, next_body,
                    ):
                        if vinf > self.config.max_encounter_vinf_mps:
                            continue
                        if next_body == target:
                            candidate = self._candidate_from_encounter(
                                state, target, vinf,
                            )
                            current = best_by_sequence.get(candidate.sequence)
                            if (
                                current is None
                                or candidate.proxy_cost_mps < current.proxy_cost_mps
                            ):
                                best_by_sequence[candidate.sequence] = candidate
                            continue
                        if len(state.sequence) + 1 >= self.config.max_bodies:
                            continue
                        body = self._body(next_body)
                        maximum_turn = maximum_flyby_turn_angle(
                            vinf,
                            body.mu_self,
                            body.safe_radius
                            + float(self.config.flyby_clearance_m),
                        )
                        for turn in _linspace(
                            -maximum_turn,
                            maximum_turn,
                            self.config.flyby_turn_samples,
                        ):
                            evidence = self._flyby_evidence(
                                state, next_body, vinf, incoming_angle, turn,
                            )
                            next_state = self._orbit_from_relative_state(
                                next_body,
                                vinf,
                                incoming_angle + turn,
                                state.sequence + (next_body,),
                                state.departure_vinf_mps,
                                max(
                                    state.max_turn_fraction,
                                    evidence.turn_fraction,
                                ),
                                state.flybys + (evidence,),
                            )
                            if next_state is not None:
                                next_frontier.append(next_state)
            if not next_frontier:
                break
            frontier = self._prune_frontier(next_frontier)

        candidates = sorted(
            best_by_sequence.values(),
            key=lambda candidate: (
                candidate.proxy_cost_mps,
                len(candidate.sequence),
                candidate.sequence,
            ),
        )
        if limit is not None:
            return candidates[:max(0, int(limit))]
        return candidates
