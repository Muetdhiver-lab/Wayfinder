# -*- coding: utf-8 -*-
"""KSP-aware translation of PyKEP MGA-1DSM genes.

The optimizer works in PyKEP's ideal patched-conic model: planets are sampled at
encounter epochs and the departure gene supplies a heliocentric excess vector.
This module keeps that model intact, then reconstructs a finite-SOI departure
state and recomputes the first DSM from that state. The result is a diagnostic
and flight-plan translation layer, not a replacement optimizer.
"""

from dataclasses import dataclass
from math import acosh, acos, atan2, pi, sinh, sqrt

import numpy as np
from numpy import cross, dot, linalg
from pykep import epoch
from pykep.core import DAY2SEC
from pykep.core import fb_vout
from scipy.linalg import norm

from _Trajectory import _lambert_problem
from _Trajectory import _lambert_v1
from _Trajectory import _lambert_v2
from _Trajectory import _max_revs
from _Trajectory import _propagate_lagrangian
from _Trajectory import _starting_soi_radius
from _Trajectory import ejection_from_gene
from _Trajectory import fast_ejection_from_gene
from _Trajectory import Kdate
from _Trajectory import Ktime
from planet_packs import PACKS


@dataclass
class ManeuverComponents:
    """Local prograde/normal/radial components for a maneuver."""

    prograde: float
    normal: float
    radial: float
    magnitude: float
    inertial: list

    def as_dict(self):
        return {
            "prograde": self.prograde,
            "normal": self.normal,
            "radial": self.radial,
            "magnitude": self.magnitude,
            "inertial": self.inertial,
        }


@dataclass
class FiniteSoiCase:
    """Finite-SOI reconstruction for one departure realization."""

    name: str
    soi_met_days: float
    soi_epoch_days: float
    soi_radius_error_m: float
    soi_target_velocity_error: float
    corrected_dsm: ManeuverComponents
    corrected_dsm_epoch_days: float
    corrected_dsm_met_days: float
    corrected_arrival_vinf: float
    corrected_total_dv_without_arrival: float
    ideal_to_corrected_dsm_delta: float

    def as_dict(self):
        return {
            "name": self.name,
            "soi_met_days": self.soi_met_days,
            "soi_epoch_days": self.soi_epoch_days,
            "soi_radius_error_m": self.soi_radius_error_m,
            "soi_target_velocity_error": self.soi_target_velocity_error,
            "corrected_dsm": self.corrected_dsm.as_dict(),
            "corrected_dsm_epoch_days": self.corrected_dsm_epoch_days,
            "corrected_dsm_met_days": self.corrected_dsm_met_days,
            "corrected_arrival_vinf": self.corrected_arrival_vinf,
            "corrected_total_dv_without_arrival": (
                self.corrected_total_dv_without_arrival
            ),
            "ideal_to_corrected_dsm_delta": self.ideal_to_corrected_dsm_delta,
        }


@dataclass
class KspTranslation:
    """KSP translation diagnostics for the first leg of a PyKEP gene."""

    arrival_epoch_days: float
    arrival_met_days: float
    dsm_epoch_days: float
    dsm_met_days: float
    ideal_dsm: ManeuverComponents
    ideal_dsm_epoch_days: float
    ideal_dsm_met_days: float
    ideal_total_dv_without_arrival: float
    selected_case: FiniteSoiCase
    split_case: FiniteSoiCase

    def as_dict(self):
        return {
            "arrival_epoch_days": self.arrival_epoch_days,
            "arrival_met_days": self.arrival_met_days,
            "dsm_epoch_days": self.dsm_epoch_days,
            "dsm_met_days": self.dsm_met_days,
            "ideal_dsm": self.ideal_dsm.as_dict(),
            "ideal_dsm_epoch_days": self.ideal_dsm_epoch_days,
            "ideal_dsm_met_days": self.ideal_dsm_met_days,
            "ideal_total_dv_without_arrival": self.ideal_total_dv_without_arrival,
            "selected_case": self.selected_case.as_dict(),
            "split_case": self.split_case.as_dict(),
        }


@dataclass
class TimelineNode:
    """One translated maneuver-node event in the mission timeline."""

    label: str
    leg_index: int
    epoch_days: float
    met_days: float
    components: ManeuverComponents
    model: str

    def as_dict(self):
        return {
            "label": self.label,
            "leg_index": self.leg_index,
            "epoch_days": self.epoch_days,
            "met_days": self.met_days,
            "components": self.components.as_dict(),
            "model": self.model,
        }


@dataclass
class TimelineFlyby:
    """One translated flyby event in the mission timeline."""

    leg_index: int
    body: str
    epoch_days: float
    met_days: float
    soi_entry_epoch_days: float
    soi_entry_met_days: float
    soi_exit_epoch_days: float
    soi_exit_met_days: float
    periapsis_radius_m: float
    periapsis_altitude_m: float
    periapsis_altitude_km: float
    radius_factor: float
    beta_angle_rad: float
    beta_angle_deg: float
    incoming_vinf: float
    outgoing_vinf: float
    hyperbola_inclination_deg: float
    flyby_plane_normal: list
    periapsis_direction: list
    model: str

    def as_dict(self):
        return {
            "leg_index": self.leg_index,
            "body": self.body,
            "epoch_days": self.epoch_days,
            "met_days": self.met_days,
            "soi_entry_epoch_days": self.soi_entry_epoch_days,
            "soi_entry_met_days": self.soi_entry_met_days,
            "soi_exit_epoch_days": self.soi_exit_epoch_days,
            "soi_exit_met_days": self.soi_exit_met_days,
            "periapsis_radius_m": self.periapsis_radius_m,
            "periapsis_altitude_m": self.periapsis_altitude_m,
            "periapsis_altitude_km": self.periapsis_altitude_km,
            "radius_factor": self.radius_factor,
            "beta_angle_rad": self.beta_angle_rad,
            "beta_angle_deg": self.beta_angle_deg,
            "incoming_vinf": self.incoming_vinf,
            "outgoing_vinf": self.outgoing_vinf,
            "hyperbola_inclination_deg": self.hyperbola_inclination_deg,
            "flyby_plane_normal": self.flyby_plane_normal,
            "periapsis_direction": self.periapsis_direction,
            "model": self.model,
        }


@dataclass
class TimelineArrival:
    """Final arrival event in the mission timeline."""

    body: str
    epoch_days: float
    met_days: float
    arrival_vinf: float
    arrival_mode: str
    model: str

    def as_dict(self):
        return {
            "body": self.body,
            "epoch_days": self.epoch_days,
            "met_days": self.met_days,
            "arrival_vinf": self.arrival_vinf,
            "arrival_mode": self.arrival_mode,
            "model": self.model,
        }


@dataclass
class KspTrajectoryTranslation:
    """Mission-level KSP translation with explicit fidelity tags."""

    first_leg: KspTranslation
    nodes: list
    flybys: list
    arrival: TimelineArrival

    def as_dict(self):
        return {
            "first_leg": self.first_leg.as_dict(),
            "nodes": [node.as_dict() for node in self.nodes],
            "flybys": [flyby.as_dict() for flyby in self.flybys],
            "arrival": self.arrival.as_dict(),
        }


def _local_components(position, velocity, delta_v):
    """Return KSP-style local components for a heliocentric maneuver."""
    prograde = np.asarray(velocity) / norm(velocity)
    plane = cross(velocity, position)
    plane = plane / norm(plane)
    radial = cross(plane, prograde)
    delta_v = np.asarray(delta_v)
    return ManeuverComponents(
        prograde=float(dot(prograde, delta_v)),
        normal=float(-dot(plane, delta_v)),
        radial=float(dot(radial, delta_v)),
        magnitude=float(norm(delta_v)),
        inertial=[float(v) for v in delta_v],
    )


def _unit_vector(vector):
    vector = np.asarray(vector, dtype=float)
    vector_norm = norm(vector)
    if vector_norm <= 1e-12:
        return np.asarray([np.nan, np.nan, np.nan], dtype=float)
    return vector / vector_norm


def _flyby_plane_geometry(incoming_vinf_vector, outgoing_vinf_vector):
    """Return pilotable flyby-plane diagnostics in the body-centered frame.

    The hyperbolic trajectory plane is spanned by v_inf_in and v_inf_out. Its
    inclination is the inclination of that hyperbola relative to the body's
    equatorial plane in the PyKEP/KSP-aligned body frame. This is the quantity
    that is easiest to compare manually in KSP once inside the flyby SOI.
    """
    incoming_hat = _unit_vector(incoming_vinf_vector)
    outgoing_hat = _unit_vector(outgoing_vinf_vector)
    plane_normal = _unit_vector(cross(incoming_hat, outgoing_hat))
    if np.any(np.isnan(plane_normal)):
        inclination = np.nan
    else:
        inclination = acos(float(np.clip(plane_normal[2], -1.0, 1.0))) * 180.0 / pi
    periapsis_velocity_hat = _unit_vector(incoming_hat + outgoing_hat)
    periapsis_direction = _unit_vector(cross(periapsis_velocity_hat, plane_normal))
    return {
        "hyperbola_inclination_deg": float(inclination),
        "flyby_plane_normal": [float(v) for v in plane_normal],
        "periapsis_direction": [float(v) for v in periapsis_direction],
    }


def _hyperbolic_soi_to_periapsis_days(body, vinf, periapsis_radius,
                                      planet_pack):
    """Approximate patched-conic time from SOI boundary to periapsis.

    This is a two-body hyperbola estimate used for timeline diagnostics. It is
    intentionally tagged as translated/approximate rather than exact KSP state
    propagation.
    """
    vinf = float(abs(vinf))
    if vinf <= 1e-9:
        return 0.0
    rsoi = _starting_soi_radius(body, planet_pack)
    rp = float(periapsis_radius)
    if rsoi <= rp:
        return 0.0
    a = float(body.mu_self) / (vinf * vinf)
    e = 1.0 + rp * vinf * vinf / float(body.mu_self)
    cosh_f = (rsoi / a + 1.0) / e
    if cosh_f < 1.0:
        return 0.0
    f_anomaly = acosh(cosh_f)
    seconds = sqrt((a ** 3) / float(body.mu_self)) * (
        e * sinh(f_anomaly) - f_anomaly
    )
    return float(seconds / DAY2SEC * PACKS[planet_pack].EDY_TO_KDY)


def _escape_to_soi(first_body, periapsis_direction, prograde_dv, normal_dv,
                   radial_dv, parking_altitude, rsoi, z_correction=0.0):
    """Propagate a local parking burn to the body's SOI in body-centered space."""
    radius = float(first_body.radius + parking_altitude)
    periapsis_direction = np.asarray(periapsis_direction, dtype=float)
    circular_hat = np.asarray([
        -periapsis_direction[1], periapsis_direction[0], 0.0
    ])
    circular_speed = sqrt(float(first_body.mu_self) / radius)
    r0 = periapsis_direction * radius
    v0 = (
        circular_speed * circular_hat
        + float(prograde_dv) * circular_hat
        + float(radial_dv) * periapsis_direction
        + np.asarray([0.0, 0.0, float(normal_dv)])
    )

    lower = 0.0
    upper = 10.0 * DAY2SEC
    for _ in range(128):
        midpoint = 0.5 * (lower + upper)
        r_mid, _ = _propagate_lagrangian(r0, v0, midpoint, first_body.mu_self)
        if norm(r_mid) < rsoi:
            lower = midpoint
        else:
            upper = midpoint
    r_soi, v_soi = _propagate_lagrangian(r0, v0, upper, first_body.mu_self)
    v_soi = np.asarray(v_soi) + np.asarray([0.0, 0.0, float(z_correction)])
    return np.asarray(r_soi), v_soi, upper


def _planar_gene_from_vinf(gene, vinf):
    planar_gene = list(gene)
    planar_gene[1] = (atan2(vinf[1], vinf[0]) % (2 * pi)) / (2 * pi)
    planar_gene[2] = 0.5
    planar_gene[3] = float(np.hypot(vinf[0], vinf[1]))
    return planar_gene


def _heliocentric_soi_state(first_body, gene, r_body_soi, v_body_soi,
                            soi_time_seconds):
    soi_epoch = float(gene[0]) + float(soi_time_seconds) / DAY2SEC
    planet_r, planet_v = first_body.eph(epoch(soi_epoch))
    return (
        np.asarray(planet_r) + np.asarray(r_body_soi),
        np.asarray(planet_v) + np.asarray(v_body_soi),
        soi_epoch,
    )


def _corrected_case(name, udp, gene, planet_pack, parking_altitude, rsoi,
                    ejection, vinf, ideal_dsm_vector, ejection_dv,
                    periapsis_direction, prograde_dv, normal_dv, radial_dv,
                    z_correction=0.0):
    first_body = udp._seq[0]
    target_body = udp._seq[1]
    edy_to_kdy = PACKS[planet_pack].EDY_TO_KDY
    times, _, _, _ = udp._decode_times_and_vinf(gene)
    first_leg_tof = float(times[0])
    dsm_met = float(gene[4]) * first_leg_tof
    arrival_epoch = float(gene[0]) + first_leg_tof

    r_soi_body, v_soi_body, soi_time = _escape_to_soi(
        first_body, periapsis_direction, prograde_dv, normal_dv, radial_dv,
        parking_altitude, rsoi, z_correction=z_correction,
    )
    r_soi_helio, v_soi_helio, soi_epoch = _heliocentric_soi_state(
        first_body, gene, r_soi_body, v_soi_body, soi_time
    )
    r_dsm, v_dsm_pre = _propagate_lagrangian(
        r_soi_helio, v_soi_helio, (float(gene[0]) + dsm_met - soi_epoch) * DAY2SEC,
        udp.common_mu,
    )
    r_target, v_target = target_body.eph(epoch(arrival_epoch))
    lambert = _lambert_problem(
        r_dsm, r_target, (first_leg_tof - dsm_met) * DAY2SEC,
        udp.common_mu, cw=False, max_revs=_max_revs(udp),
    )
    corrected_departure = np.asarray(_lambert_v1(lambert))
    corrected_arrival = np.asarray(_lambert_v2(lambert))
    corrected_dsm_vector = corrected_departure - np.asarray(v_dsm_pre)
    corrected_dsm = _local_components(r_dsm, v_dsm_pre, corrected_dsm_vector)
    target_vinf = corrected_arrival - np.asarray(v_target)

    target_velocity_error = norm(np.asarray(v_soi_body) - np.asarray(vinf))
    return FiniteSoiCase(
        name=name,
        soi_met_days=float(soi_time / DAY2SEC * edy_to_kdy),
        soi_epoch_days=float(soi_epoch * edy_to_kdy),
        soi_radius_error_m=float(abs(norm(r_soi_body) - rsoi)),
        soi_target_velocity_error=float(target_velocity_error),
        corrected_dsm=corrected_dsm,
        corrected_dsm_epoch_days=float((float(gene[0]) + dsm_met) * edy_to_kdy),
        corrected_dsm_met_days=float(dsm_met * edy_to_kdy),
        corrected_arrival_vinf=float(norm(target_vinf)),
        corrected_total_dv_without_arrival=float(
            ejection_dv + corrected_dsm.magnitude
        ),
        ideal_to_corrected_dsm_delta=float(
            norm(corrected_dsm_vector - np.asarray(ideal_dsm_vector))
        ),
    )


def translate_first_leg(udp, gene, planet_pack="Vanilla",
                        parking_altitude=100000.0):
    """Translate the first leg of an MGA-1DSM gene into KSP-aware diagnostics."""
    if udp.n_legs < 1:
        raise ValueError("Expected at least one leg")
    edy_to_kdy = PACKS[planet_pack].EDY_TO_KDY
    times, vinfx, vinfy, vinfz = udp._decode_times_and_vinf(gene)
    vinf = np.asarray([vinfx, vinfy, vinfz], dtype=float)
    rsoi = _starting_soi_radius(udp._seq[0], planet_pack)
    ejection = ejection_from_gene(udp, gene, rsoi, parking_altitude)

    r0, v_planet0 = udp._seq[0].eph(epoch(gene[0]))
    v0 = np.asarray(v_planet0) + vinf
    r_ideal_dsm, v_ideal_dsm_pre = _propagate_lagrangian(
        r0, v0, float(gene[4]) * float(times[0]) * DAY2SEC, udp.common_mu
    )
    r_target, _ = udp._seq[1].eph(epoch(float(gene[0]) + float(times[0])))
    lambert = _lambert_problem(
        r_ideal_dsm, r_target,
        (1.0 - float(gene[4])) * float(times[0]) * DAY2SEC,
        udp.common_mu, cw=False, max_revs=_max_revs(udp),
    )
    ideal_dsm_vector = np.asarray(_lambert_v1(lambert)) - np.asarray(v_ideal_dsm_pre)
    ideal_dsm = _local_components(r_ideal_dsm, v_ideal_dsm_pre, ideal_dsm_vector)

    selected_case = _corrected_case(
        "selected",
        udp, gene, planet_pack, parking_altitude, rsoi,
        ejection, vinf, ideal_dsm_vector, ejection["dv"],
        ejection["periapsis_direction"],
        ejection["parking_node_prograde_dv"],
        ejection["parking_node_normal_dv"],
        ejection["parking_node_radial_dv"],
    )

    planar_gene = _planar_gene_from_vinf(gene, vinf)
    planar = fast_ejection_from_gene(
        udp._seq[0], planar_gene, parking_altitude, rsoi
    )
    split_case = _corrected_case(
        "split_soi",
        udp, gene, planet_pack, parking_altitude, rsoi,
        ejection, vinf, ideal_dsm_vector, planar["parking_node_dv"] + abs(vinf[2]),
        planar["periapsis_direction"],
        planar["parking_node_prograde_dv"],
        planar["parking_node_normal_dv"],
        planar["parking_node_radial_dv"],
        z_correction=vinf[2],
    )

    return KspTranslation(
        arrival_epoch_days=float((float(gene[0]) + float(times[0])) * edy_to_kdy),
        arrival_met_days=float(float(times[0]) * edy_to_kdy),
        dsm_epoch_days=float(
            (float(gene[0]) + float(gene[4]) * float(times[0])) * edy_to_kdy
        ),
        dsm_met_days=float(float(gene[4]) * float(times[0]) * edy_to_kdy),
        ideal_dsm=ideal_dsm,
        ideal_dsm_epoch_days=float(
            (float(gene[0]) + float(gene[4]) * float(times[0])) * edy_to_kdy
        ),
        ideal_dsm_met_days=float(float(gene[4]) * float(times[0]) * edy_to_kdy),
        ideal_total_dv_without_arrival=float(ejection["dv"] + ideal_dsm.magnitude),
        selected_case=selected_case,
        split_case=split_case,
    )


def translate_mga_trajectory(
    udp,
    gene,
    planet_pack="Vanilla",
    parking_altitude=100000.0,
    arrival_mode="flyby",
):
    """Translate an MGA-1DSM gene into a mission-level KSP timeline.

    The first leg receives the finite-SOI correction implemented by
    :func:`translate_first_leg`. Later flybys are converted into KSP timeline
    diagnostics with SOI-entry/Pe/SOI-exit epochs and v-infinity values. The
    post-flyby DSM vectors currently remain PyKEP-ideal fallback events and are
    explicitly tagged as such.
    """
    first_leg = translate_first_leg(
        udp, gene, planet_pack=planet_pack, parking_altitude=parking_altitude
    )
    edy_to_kdy = PACKS[planet_pack].EDY_TO_KDY
    times, vinfx, vinfy, vinfz = udp._decode_times_and_vinf(gene)
    vinf = np.asarray([vinfx, vinfy, vinfz], dtype=float)
    departure_epoch_days = float(gene[0]) * edy_to_kdy

    t_p = []
    r_p = []
    v_p = []
    for i, planet in enumerate(udp._seq):
        t_i = epoch(float(gene[0] + sum(times[0:i])))
        r_i, v_i = planet.eph(t_i)
        t_p.append(t_i)
        r_p.append(np.asarray(r_i, dtype=float))
        v_p.append(np.asarray(v_i, dtype=float))

    nodes = []
    flybys = []

    v0 = v_p[0] + vinf
    dsm_fraction = float(gene[4])
    r, v = _propagate_lagrangian(
        r_p[0], v0, dsm_fraction * times[0] * DAY2SEC, udp.common_mu
    )
    lambert = _lambert_problem(
        r,
        r_p[1],
        (1.0 - dsm_fraction) * times[0] * DAY2SEC,
        udp.common_mu,
        cw=False,
        max_revs=_max_revs(udp),
    )
    v_end_l = np.asarray(_lambert_v2(lambert), dtype=float)

    for i in range(1, udp.n_legs):
        body = udp._seq[i]
        radius_factor = float(gene[7 + (i - 1) * 4])
        periapsis_radius = radius_factor * float(body.radius)
        periapsis_altitude = periapsis_radius - float(body.radius)
        beta = float(gene[6 + (i - 1) * 4])
        flyby_epoch_days = float(t_p[i].mjd2000 * edy_to_kdy)
        incoming_vinf_vector = v_end_l - v_p[i]
        incoming_vinf = float(norm(incoming_vinf_vector))
        v_out = np.asarray(
            fb_vout(v_end_l, v_p[i], periapsis_radius, beta, body.mu_self),
            dtype=float,
        )
        outgoing_vinf = float(norm(v_out - v_p[i]))
        flyby_geometry = _flyby_plane_geometry(
            incoming_vinf_vector, v_out - v_p[i]
        )
        entry_offset = _hyperbolic_soi_to_periapsis_days(
            body, incoming_vinf, periapsis_radius, planet_pack
        )
        exit_offset = _hyperbolic_soi_to_periapsis_days(
            body, outgoing_vinf, periapsis_radius, planet_pack
        )
        flybys.append(
            TimelineFlyby(
                leg_index=i,
                body=body.name,
                epoch_days=flyby_epoch_days,
                met_days=flyby_epoch_days - departure_epoch_days,
                soi_entry_epoch_days=flyby_epoch_days - entry_offset,
                soi_entry_met_days=flyby_epoch_days - entry_offset
                - departure_epoch_days,
                soi_exit_epoch_days=flyby_epoch_days + exit_offset,
                soi_exit_met_days=flyby_epoch_days + exit_offset
                - departure_epoch_days,
                periapsis_radius_m=float(periapsis_radius),
                periapsis_altitude_m=float(periapsis_altitude),
                periapsis_altitude_km=float(periapsis_altitude / 1000.0),
                radius_factor=radius_factor,
                beta_angle_rad=beta,
                beta_angle_deg=float(beta * 180.0 / pi),
                incoming_vinf=incoming_vinf,
                outgoing_vinf=outgoing_vinf,
                hyperbola_inclination_deg=flyby_geometry[
                    "hyperbola_inclination_deg"
                ],
                flyby_plane_normal=flyby_geometry["flyby_plane_normal"],
                periapsis_direction=flyby_geometry["periapsis_direction"],
                model="ksp_flyby_translated",
            )
        )

        leg_dsm_fraction = float(gene[8 + (i - 1) * 4])
        r, v = _propagate_lagrangian(
            r_p[i],
            v_out,
            leg_dsm_fraction * times[i] * DAY2SEC,
            udp.common_mu,
        )
        lambert = _lambert_problem(
            r,
            r_p[i + 1],
            (1.0 - leg_dsm_fraction) * times[i] * DAY2SEC,
            udp.common_mu,
            cw=False,
            max_revs=_max_revs(udp),
        )
        v_beg_l = np.asarray(_lambert_v1(lambert), dtype=float)
        v_end_l = np.asarray(_lambert_v2(lambert), dtype=float)
        dsm_vector = v_beg_l - np.asarray(v, dtype=float)
        dsm_epoch_days = float(
            (t_p[i].mjd2000 + leg_dsm_fraction * times[i]) * edy_to_kdy
        )
        nodes.append(
            TimelineNode(
                label=f"dsm_leg_{i + 1}_{body.name}_to_{udp._seq[i + 1].name}",
                leg_index=i + 1,
                epoch_days=dsm_epoch_days,
                met_days=dsm_epoch_days - departure_epoch_days,
                components=_local_components(r, v, dsm_vector),
                model="pykep_ideal_fallback",
            )
        )

    arrival_epoch_days = float(t_p[-1].mjd2000 * edy_to_kdy)
    arrival = TimelineArrival(
        body=udp._seq[-1].name,
        epoch_days=arrival_epoch_days,
        met_days=arrival_epoch_days - departure_epoch_days,
        arrival_vinf=float(norm(v_end_l - v_p[-1])),
        arrival_mode=arrival_mode,
        model="pykep_arrival_fallback",
    )
    return KspTrajectoryTranslation(
        first_leg=first_leg,
        nodes=nodes,
        flybys=flybys,
        arrival=arrival,
    )


def format_translation_summary(translation, planet_pack="Vanilla"):
    """Return human-readable lines for transx-style output."""
    lines = []
    lines.append("KSP translation diagnostics:")
    lines.append(
        "  translated TOF:       "
        + f"T+{translation.arrival_met_days:.4f} d, "
        + f"arrival {translation.arrival_epoch_days:.4f} KUT "
        + f"({Kdate(translation.arrival_epoch_days, planet_pack)})"
    )
    lines.append(
        "  translated DSM epoch: "
        + f"T+{translation.dsm_met_days:.4f} d, "
        + f"{translation.dsm_epoch_days:.4f} KUT "
        + f"({Kdate(translation.dsm_epoch_days, planet_pack)})"
    )
    lines.append(
        "  ideal DSM:            "
        + f"{translation.ideal_dsm.prograde:.3f} m/s prograde, "
        + f"{translation.ideal_dsm.normal:.3f} m/s normal, "
        + f"{translation.ideal_dsm.radial:.3f} m/s radial "
        + f"(|dV|={translation.ideal_dsm.magnitude:.3f} m/s)"
    )
    for case in (translation.selected_case, translation.split_case):
        lines.append(
            f"  {case.name} SOI:          "
            + f"T+{case.soi_met_days:.4f} d, "
            + f"{case.soi_epoch_days:.4f} KUT "
            + f"({Kdate(case.soi_epoch_days, planet_pack)})"
        )
        lines.append(
            f"  {case.name} corrected DSM:"
            + f" {case.corrected_dsm.prograde:.3f} P,"
            + f" {case.corrected_dsm.normal:.3f} N,"
            + f" {case.corrected_dsm.radial:.3f} R"
            + f" (|dV|={case.corrected_dsm.magnitude:.3f} m/s)"
        )
        lines.append(
            f"    operational DV w/o arrival: "
            + f"{case.corrected_total_dv_without_arrival:.3f} m/s"
            + f" (+{case.ideal_to_corrected_dsm_delta:.3f} m/s DSM shift)"
        )
    return lines
