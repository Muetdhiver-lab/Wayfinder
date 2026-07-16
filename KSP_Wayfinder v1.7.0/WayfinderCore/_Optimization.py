# -*- coding: utf-8 -*-
"""Optimization decorators and telemetry helpers for Wayfinder."""

import ast
import logging
import re
from math import acos, cos, log, pi, sin, sqrt

from _Trajectory import fast_ejection_from_gene


logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


def alpha_leg_tofs(dv, n_legs):
    """Decode pykep mga_1dsm alpha-encoded leg times of flight."""
    weights = [-log(dv[5 + 4 * i]) for i in range(n_legs)]
    alpha_sum = sum(weights)
    return [dv[-1] * weight / alpha_sum for weight in weights]


def direct_leg_tofs(dv, n_legs):
    """Decode pykep mga_1dsm direct-encoded leg times of flight."""
    return [dv[5 + 4 * i] for i in range(n_legs)]


def leg_tofs(dv, n_legs, tof_encoding):
    """Decode leg times of flight for supported pykep mga_1dsm encodings."""
    if tof_encoding == "alpha":
        return alpha_leg_tofs(dv, n_legs)
    if tof_encoding == "direct":
        return direct_leg_tofs(dv, n_legs)
    raise ValueError("Unsupported tof_encoding: " + str(tof_encoding))


def alpha_gene_to_direct_gene(alpha_gene, n_legs):
    """Convert an alpha-encoded mga_1dsm gene to direct encoding."""
    direct_gene = list(alpha_gene[:5])
    decoded_tofs = alpha_leg_tofs(alpha_gene, n_legs)
    for leg_index in range(n_legs):
        direct_gene.append(decoded_tofs[leg_index])
        if leg_index < n_legs - 1:
            offset = 6 + 4 * leg_index
            direct_gene.extend(alpha_gene[offset:offset + 3])
    return direct_gene


def direct_tof_bounds_from_leg_tofs(decoded_tofs, wiggle):
    """Build pykep direct-mode per-leg TOF bounds around decoded leg TOFs."""
    if wiggle < 0:
        raise ValueError("wiggle must be non-negative")
    bounds = []
    for tof in decoded_tofs:
        lower = max(1e-6, float(tof) - float(wiggle))
        upper = max(lower + 1e-6, float(tof) + float(wiggle))
        bounds.append([lower, upper])
    return bounds


def _sequence_from_problem(problem):
    sequence_rgx = re.compile(r"\[.*?\]")
    bracket = sequence_rgx.search(problem.inner_problem.get_extra_info())
    return ast.literal_eval(bracket.group())


class WayfinderFitnessDecorator:
    """Decorate pykep's mga_1dsm fitness with Wayfinder-specific penalties."""

    def __init__(
        self, planet_pack, bodies_by_name, soi_radius_by_name,
        ejection_altitude, tof_encoding="alpha", ejection_model="approximate",
    ):
        self.planet_pack = planet_pack
        self.bodies_by_name = bodies_by_name
        self.soi_radius_by_name = soi_radius_by_name
        self.ejection_altitude = ejection_altitude
        self.tof_encoding = tof_encoding
        self.ejection_model = ejection_model

    def __call__(self, orig_fitness_function):
        def new_fitness_function(problem, dv):
            fitness = self._fitness_with_ejection_cost(problem, dv, orig_fitness_function)
            sequence = _sequence_from_problem(problem)
            leg_tofs_values = leg_tofs(dv, len(sequence) - 1, self.tof_encoding)
            if self.ejection_model != "vector_3d":
                fitness += self._inclination_penalty(dv)
            fitness += self._legacy_kerbin_resonance_penalty(sequence, leg_tofs_values)
            return fitness

        return new_fitness_function

    def _fitness_with_ejection_cost(self, problem, dv, orig_fitness_function):
        sequence = _sequence_from_problem(problem)
        first_body = self.bodies_by_name[sequence[0]]

        theta = 2 * pi * dv[1]
        phi = acos(2 * dv[2] - 1) - pi / 2
        vinfx = dv[3] * cos(phi) * cos(theta)
        vinfy = dv[3] * cos(phi) * sin(theta)
        vinfz = dv[3] * sin(phi)
        vinfxy = sqrt(vinfy**2 + vinfx**2)

        if self.ejection_model == "vector_3d":
            r_soi = self.soi_radius_by_name[first_body.name]
            ejection_dv = fast_ejection_from_gene(
                first_body, dv, self.ejection_altitude, r_soi
            )["dv"]
            return orig_fitness_function(problem, dv) + ejection_dv - dv[3]
        if self.ejection_model != "approximate":
            raise ValueError("Unsupported ejection model: " + str(self.ejection_model))

        r_soi = self.soi_radius_by_name[first_body.name]
        v0 = sqrt(first_body.mu_self / (first_body.radius + self.ejection_altitude))
        vxy_ob = (
            sqrt(
                vinfxy**2
                + 2 * first_body.mu_self / (first_body.radius + self.ejection_altitude)
                - 2 * first_body.mu_self / r_soi
            )
            - v0
        )
        v_ob = sqrt(vxy_ob**2 + vinfz**2)

        if v_ob < (vxy_ob + abs(vinfz)):
            return orig_fitness_function(problem, dv) + v_ob - dv[3]
        return orig_fitness_function(problem, dv) + vxy_ob + abs(vinfz) - dv[3]

    def _inclination_penalty(self, dv):
        theta = 2 * pi * dv[1]
        phi = acos(2 * dv[2] - 1) - pi / 2
        vinfx = dv[3] * cos(phi) * cos(theta)
        vinfy = dv[3] * cos(phi) * sin(theta)
        vinfz = dv[3] * sin(phi)
        vinfxy = sqrt(vinfy**2 + vinfx**2)
        if vinfxy == 0:
            return 0 if vinfz == 0 else 1e12
        return max(0, abs(vinfz / vinfxy) - 0.25) * 10000

    def _legacy_kerbin_resonance_penalty(self, sequence, leg_tofs):
        if self.planet_pack == "Vanilla":
            return self._vanilla_kerbin_resonance_penalty(sequence, leg_tofs)
        if self.planet_pack == "JNSQ":
            return self._jnsq_kerbin_resonance_penalty(sequence, leg_tofs)
        logger.warning("No resonance penalty configured for planet pack %s", self.planet_pack)
        return 0

    def _vanilla_kerbin_resonance_penalty(self, sequence, leg_tofs):
        penalty = 0
        kerbin_year = 426 / 4.0
        half_kerbin_year_window = kerbin_year * 0.52
        crop_window = 5 / 4.0

        if sequence[0] == "Kerbin" and sequence[0] == sequence[1]:
            penalty += 10000 * (
                max(half_kerbin_year_window, abs(leg_tofs[0] - kerbin_year))
                - half_kerbin_year_window
            )

        for planet, next_planet, tof in zip(sequence[1:], sequence[2:], leg_tofs[1:]):
            if planet == "Kerbin" and next_planet == planet:
                penalty += 10000 * (max(crop_window, abs(tof - kerbin_year * 2)) - crop_window)
        return penalty

    def _jnsq_kerbin_resonance_penalty(self, sequence, leg_tofs):
        penalty = 0
        kerbin_two_years = 2 * 365 / 2
        kerbin_window = 5 / 2.0
        for planet, next_planet, tof in zip(sequence, sequence[1:], leg_tofs):
            if planet == "Kerbin" and next_planet == planet:
                penalty += 10000 * (max(kerbin_window, abs(tof - kerbin_two_years)) - kerbin_window)
        return penalty
