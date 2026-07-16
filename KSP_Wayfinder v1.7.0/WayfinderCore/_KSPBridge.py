# -*- coding: utf-8 -*-
"""Optional kRPC bridge for injecting Wayfinder flight plans into KSP.

This module is intentionally optional: importing Wayfinder should not require
the ``krpc`` Python package or a running game. The bridge consumes the
KSP-aware translation layer and turns it into concrete maneuver nodes.

Scope note
----------
This bridge is deliberately an injector, not an autopilot or a flyby
corrector. Multi-flyby KSP patched conics are fragile over long arcs; Wayfinder
provides the node sequence plus pilot-facing flyby targets such as Pe altitude
and hyperbola inclination, while leg-by-leg correction remains manual.
"""

from dataclasses import dataclass
from typing import List, Optional

from planet_packs import PACKS

from _KSPTranslator import translate_mga_trajectory
from _Trajectory import _starting_soi_radius
from _Trajectory import ejection_from_gene
from _Trajectory import fast_ejection_from_gene
from _Trajectory import Kdate

import math
import numpy as np
from pykep.core import DAY2SEC


@dataclass
class KspNodePlan:
    """One maneuver node in KSP UT seconds and local node components."""

    label: str
    epoch_days: float
    met_days: float
    ut_seconds: float
    prograde: float
    normal: float
    radial: float
    model: str = "ksp"
    leg_index: Optional[int] = None

    @property
    def magnitude(self):
        return float(
            (self.prograde ** 2 + self.normal ** 2 + self.radial ** 2) ** 0.5
        )

    def as_dict(self):
        return {
            "label": self.label,
            "epoch_days": self.epoch_days,
            "met_days": self.met_days,
            "ut_seconds": self.ut_seconds,
            "prograde": self.prograde,
            "normal": self.normal,
            "radial": self.radial,
            "magnitude": self.magnitude,
            "model": self.model,
            "leg_index": self.leg_index,
        }


@dataclass
class KspFlybyPlan:
    """One planned gravity-assist event."""

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
class KspArrivalPlan:
    """Final target-arrival event."""

    body: str
    epoch_days: float
    met_days: float
    arrival_vinf: float
    arrival_mode: str = "flyby"
    model: str = "pykep_arrival_fallback"

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
class KspFlightPlan:
    """A KSP-actionable maneuver-node plan plus mission timeline."""

    planet_pack: str
    strategy: str
    departure_epoch_days: float
    arrival_epoch_days: float
    arrival_met_days: float
    nodes: List[KspNodePlan]
    flybys: List[KspFlybyPlan]
    arrival: Optional[KspArrivalPlan]
    ideal_total_dv_without_arrival: float
    operational_total_dv_without_arrival: float

    def as_dict(self):
        return {
            "planet_pack": self.planet_pack,
            "strategy": self.strategy,
            "departure_epoch_days": self.departure_epoch_days,
            "arrival_epoch_days": self.arrival_epoch_days,
            "arrival_met_days": self.arrival_met_days,
            "nodes": [node.as_dict() for node in self.nodes],
            "flybys": [flyby.as_dict() for flyby in self.flybys],
            "arrival": None if self.arrival is None else self.arrival.as_dict(),
            "ideal_total_dv_without_arrival": self.ideal_total_dv_without_arrival,
            "operational_total_dv_without_arrival": (
                self.operational_total_dv_without_arrival
            ),
        }

    def format(self, title=None, metadata=None):
        """Return a pilot-facing text flight plan with node DV breakdown."""
        return format_ksp_flight_plan(self, title=title, metadata=metadata)


def _kday_to_seconds(planet_pack):
    return DAY2SEC / PACKS[planet_pack].EDY_TO_KDY


def _fmt_days(value):
    return f"{float(value):.4f}"


def _fmt_dv(value):
    return f"{float(value):.3f}"


def _fmt_vector(values, digits=5):
    return "[" + ", ".join(f"{float(v):.{digits}f}" for v in values) + "]"


def format_ksp_flight_plan(plan, title=None, metadata=None):
    """Format a KSP-actionable flight plan for manual flying.

    The output intentionally distinguishes the injected/operational maneuver
    node budget from arrival v-infinity. For ``arrival_mode='flyby'`` the final
    v-infinity is a diagnostic, not a burn.
    """
    metadata = metadata or {}
    lines = []
    lines.append(title or "KSP flight plan")
    lines.append("=" * len(lines[-1]))
    if metadata:
        for key in (
            "run_id",
            "sequence",
            "sequence_short_name",
            "objective_dv",
            "runtime_seconds",
            "result_t0",
            "result_tof",
        ):
            if key in metadata and metadata[key] is not None:
                label = key.replace("_", " ")
                value = metadata[key]
                if key in ("objective_dv",):
                    value = f"{float(value):.3f} m/s"
                elif key in ("runtime_seconds",):
                    value = f"{float(value):.3f} s"
                elif key in ("result_t0", "result_tof"):
                    value = f"{float(value):.4f} days"
                lines.append(f"{label:<24}: {value}")
    lines.append(f"{'planet pack':<24}: {plan.planet_pack}")
    lines.append(f"{'strategy':<24}: {plan.strategy}")
    lines.append(
        f"{'departure':<24}: {_fmt_days(plan.departure_epoch_days)} KUT "
        f"({Kdate(plan.departure_epoch_days, plan.planet_pack)})"
    )
    lines.append(
        f"{'arrival':<24}: {_fmt_days(plan.arrival_epoch_days)} KUT "
        f"({Kdate(plan.arrival_epoch_days, plan.planet_pack)}), "
        f"T+{_fmt_days(plan.arrival_met_days)} d"
    )
    lines.append(
        f"{'first-leg ideal DV':<24}: "
        f"{_fmt_dv(plan.ideal_total_dv_without_arrival)} m/s"
    )
    lines.append(
        f"{'operational node DV':<24}: "
        f"{_fmt_dv(plan.operational_total_dv_without_arrival)} m/s"
    )

    by_model = {}
    by_leg = {}
    for node in plan.nodes:
        by_model[node.model] = by_model.get(node.model, 0.0) + node.magnitude
        leg_key = "unknown" if node.leg_index is None else f"leg {node.leg_index}"
        by_leg[leg_key] = by_leg.get(leg_key, 0.0) + node.magnitude

    lines.append("")
    lines.append("DV breakdown")
    lines.append("------------")
    for model, dv in sorted(by_model.items()):
        lines.append(f"  by model {model:<28}: {_fmt_dv(dv)} m/s")
    for leg, dv in sorted(by_leg.items()):
        lines.append(f"  by {leg:<34}: {_fmt_dv(dv)} m/s")

    lines.append("")
    lines.append("Maneuver nodes")
    lines.append("--------------")
    for index, node in enumerate(plan.nodes, start=1):
        lines.append(
            f"{index:02d}. {node.label} [{node.model}]"
        )
        lines.append(
            f"    KUT {_fmt_days(node.epoch_days)} "
            f"({Kdate(node.epoch_days, plan.planet_pack)}), "
            f"T+{_fmt_days(node.met_days)} d"
        )
        lines.append(
            f"    P {node.prograde:+.3f} m/s, "
            f"N {node.normal:+.3f} m/s, "
            f"R {node.radial:+.3f} m/s, "
            f"|dV| {_fmt_dv(node.magnitude)} m/s"
        )

    if plan.flybys:
        lines.append("")
        lines.append("Flyby targets")
        lines.append("-------------")
        for flyby in plan.flybys:
            lines.append(
                f"- {flyby.body} flyby, leg {flyby.leg_index} "
                f"[{flyby.model}]"
            )
            lines.append(
                f"    SOI entry {_fmt_days(flyby.soi_entry_epoch_days)} KUT "
                f"({Kdate(flyby.soi_entry_epoch_days, plan.planet_pack)}), "
                f"T+{_fmt_days(flyby.soi_entry_met_days)} d"
            )
            lines.append(
                f"    Pe        {_fmt_days(flyby.epoch_days)} KUT "
                f"({Kdate(flyby.epoch_days, plan.planet_pack)}), "
                f"T+{_fmt_days(flyby.met_days)} d"
            )
            lines.append(
                f"    SOI exit  {_fmt_days(flyby.soi_exit_epoch_days)} KUT "
                f"({Kdate(flyby.soi_exit_epoch_days, plan.planet_pack)}), "
                f"T+{_fmt_days(flyby.soi_exit_met_days)} d"
            )
            lines.append(
                f"    Pe altitude {flyby.periapsis_altitude_km:.3f} km, "
                f"radius {flyby.periapsis_radius_m / 1000.0:.3f} km"
            )
            lines.append(
                f"    hyperbola inclination {flyby.hyperbola_inclination_deg:.3f} deg, "
                f"beta {flyby.beta_angle_deg:.3f} deg"
            )
            lines.append(
                f"    incoming v_inf {_fmt_dv(flyby.incoming_vinf)} m/s, "
                f"outgoing v_inf {_fmt_dv(flyby.outgoing_vinf)} m/s"
            )
            lines.append(
                f"    plane normal {_fmt_vector(flyby.flyby_plane_normal)}, "
                f"periapsis dir {_fmt_vector(flyby.periapsis_direction)}"
            )

    if plan.arrival is not None:
        lines.append("")
        lines.append("Arrival")
        lines.append("-------")
        lines.append(
            f"{plan.arrival.body}: {_fmt_days(plan.arrival.epoch_days)} KUT "
            f"({Kdate(plan.arrival.epoch_days, plan.planet_pack)}), "
            f"T+{_fmt_days(plan.arrival.met_days)} d"
        )
        lines.append(
            f"arrival mode {plan.arrival.arrival_mode}, "
            f"arrival v_inf {_fmt_dv(plan.arrival.arrival_vinf)} m/s"
        )
        if plan.arrival.arrival_mode == "flyby":
            lines.append("arrival v_inf is diagnostic only; no capture burn is included.")

    return "\n".join(lines)


def _epoch_days_to_ut(epoch_days, planet_pack):
    return float(epoch_days) * _kday_to_seconds(planet_pack)


def _planar_gene_from_vinf(gene, vinf):
    planar_gene = list(gene)
    planar_gene[1] = (math.atan2(vinf[1], vinf[0]) % (2 * math.pi)) / (
        2 * math.pi
    )
    planar_gene[2] = 0.5
    planar_gene[3] = float(np.hypot(vinf[0], vinf[1]))
    return planar_gene


def build_first_leg_flight_plan(
    udp,
    gene,
    planet_pack="Vanilla",
    parking_altitude=100000.0,
    strategy="selected",
    corrected_dsm=True,
    arrival_mode="flyby",
):
    """Build a maneuver-node plan for the first leg.

    Parameters
    ----------
    strategy:
        ``"selected"`` uses the selected/direct ejection from the PyKEP gene.
        ``"split_soi"`` uses a planar parking burn, a normal SOI correction,
        then the split finite-SOI corrected DSM.
    corrected_dsm:
        If true, inject the finite-SOI corrected DSM. If false, inject the
        ideal PyKEP DSM for comparison/debugging.
    """
    if strategy not in ("selected", "split_soi"):
        raise ValueError("strategy must be 'selected' or 'split_soi'")

    times, vinfx, vinfy, vinfz = udp._decode_times_and_vinf(gene)
    vinf = np.asarray([vinfx, vinfy, vinfz], dtype=float)
    rsoi = _starting_soi_radius(udp._seq[0], planet_pack)
    ejection = ejection_from_gene(udp, gene, rsoi, parking_altitude)
    mission_translation = translate_mga_trajectory(
        udp,
        gene,
        planet_pack=planet_pack,
        parking_altitude=parking_altitude,
        arrival_mode=arrival_mode,
    )
    translation = mission_translation.first_leg
    departure_epoch_days = float(gene[0]) * PACKS[planet_pack].EDY_TO_KDY
    nodes = []

    if strategy == "selected":
        nodes.append(
            KspNodePlan(
                label="departure",
                epoch_days=departure_epoch_days,
                met_days=0.0,
                ut_seconds=_epoch_days_to_ut(departure_epoch_days, planet_pack),
                prograde=float(ejection["parking_node_prograde_dv"]),
                normal=float(ejection["parking_node_normal_dv"]),
                radial=float(ejection["parking_node_radial_dv"]),
                model="ksp_finite_soi",
                leg_index=1,
            )
        )
        case = translation.selected_case
    else:
        planar_gene = _planar_gene_from_vinf(gene, vinf)
        planar = fast_ejection_from_gene(
            udp._seq[0], planar_gene, parking_altitude, rsoi
        )
        nodes.append(
            KspNodePlan(
                label="departure_planar",
                epoch_days=departure_epoch_days,
                met_days=0.0,
                ut_seconds=_epoch_days_to_ut(departure_epoch_days, planet_pack),
                prograde=float(planar["parking_node_prograde_dv"]),
                normal=float(planar["parking_node_normal_dv"]),
                radial=float(planar["parking_node_radial_dv"]),
                model="ksp_finite_soi",
                leg_index=1,
            )
        )
        case = translation.split_case
        nodes.append(
            KspNodePlan(
                label="soi_normal_correction",
                epoch_days=float(case.soi_epoch_days),
                met_days=float(case.soi_met_days),
                ut_seconds=_epoch_days_to_ut(case.soi_epoch_days, planet_pack),
                prograde=0.0,
                normal=float(vinf[2]),
                radial=0.0,
                model="ksp_finite_soi",
                leg_index=1,
            )
        )

    if corrected_dsm:
        dsm = case.corrected_dsm
        dsm_label = "dsm_corrected"
    else:
        dsm = translation.ideal_dsm
        dsm_label = "dsm_ideal"
    nodes.append(
        KspNodePlan(
            label=dsm_label,
            epoch_days=float(translation.dsm_epoch_days),
            met_days=float(translation.dsm_met_days),
            ut_seconds=_epoch_days_to_ut(translation.dsm_epoch_days, planet_pack),
            prograde=float(dsm.prograde),
            normal=float(dsm.normal),
            radial=float(dsm.radial),
            model=(
                "ksp_finite_soi_corrected"
                if corrected_dsm
                else "pykep_ideal"
            ),
            leg_index=1,
        )
    )

    flybys = [_bridge_flyby_from_translation(flyby) for flyby in mission_translation.flybys]
    post_flyby_nodes = [
        _bridge_node_from_translation(node, planet_pack)
        for node in mission_translation.nodes
    ]
    arrival = _bridge_arrival_from_translation(mission_translation.arrival)
    nodes.extend(post_flyby_nodes)
    operational_total_dv = float(sum(node.magnitude for node in nodes))

    return KspFlightPlan(
        planet_pack=planet_pack,
        strategy=strategy,
        departure_epoch_days=departure_epoch_days,
        arrival_epoch_days=float(arrival.epoch_days),
        arrival_met_days=float(arrival.met_days),
        nodes=nodes,
        flybys=flybys,
        arrival=arrival,
        ideal_total_dv_without_arrival=float(
            translation.ideal_total_dv_without_arrival
        ),
        operational_total_dv_without_arrival=operational_total_dv,
    )


def build_ksp_flight_plan(*args, **kwargs):
    """Build a mission-level KSP flight plan.

    This is the preferred public entry point. ``build_first_leg_flight_plan``
    is retained as a compatibility alias for older tests and scripts.
    """
    return build_first_leg_flight_plan(*args, **kwargs)


def _bridge_node_from_translation(node, planet_pack):
    components = node.components
    return KspNodePlan(
        label=node.label,
        epoch_days=float(node.epoch_days),
        met_days=float(node.met_days),
        ut_seconds=_epoch_days_to_ut(node.epoch_days, planet_pack),
        prograde=float(components.prograde),
        normal=float(components.normal),
        radial=float(components.radial),
        model=node.model,
        leg_index=int(node.leg_index),
    )


def _bridge_flyby_from_translation(flyby):
    return KspFlybyPlan(
        leg_index=int(flyby.leg_index),
        body=flyby.body,
        epoch_days=float(flyby.epoch_days),
        met_days=float(flyby.met_days),
        soi_entry_epoch_days=float(flyby.soi_entry_epoch_days),
        soi_entry_met_days=float(flyby.soi_entry_met_days),
        soi_exit_epoch_days=float(flyby.soi_exit_epoch_days),
        soi_exit_met_days=float(flyby.soi_exit_met_days),
        periapsis_radius_m=float(flyby.periapsis_radius_m),
        periapsis_altitude_m=float(flyby.periapsis_altitude_m),
        periapsis_altitude_km=float(flyby.periapsis_altitude_km),
        radius_factor=float(flyby.radius_factor),
        beta_angle_rad=float(flyby.beta_angle_rad),
        beta_angle_deg=float(flyby.beta_angle_deg),
        incoming_vinf=float(flyby.incoming_vinf),
        outgoing_vinf=float(flyby.outgoing_vinf),
        hyperbola_inclination_deg=float(flyby.hyperbola_inclination_deg),
        flyby_plane_normal=list(flyby.flyby_plane_normal),
        periapsis_direction=list(flyby.periapsis_direction),
        model=flyby.model,
    )


def _bridge_arrival_from_translation(arrival):
    return KspArrivalPlan(
        body=arrival.body,
        epoch_days=float(arrival.epoch_days),
        met_days=float(arrival.met_days),
        arrival_vinf=float(arrival.arrival_vinf),
        arrival_mode=arrival.arrival_mode,
        model=arrival.model,
    )


class KspKrpcClient:
    """Small kRPC client wrapper for maneuver-node injection only."""

    def __init__(self, name="Wayfinder", address="127.0.0.1",
                 rpc_port=50000, stream_port=None):
        self.name = name
        self.address = address
        self.rpc_port = rpc_port
        self.stream_port = stream_port
        self.connection = None

    def connect(self):
        import krpc

        self.connection = krpc.connect(
            name=self.name,
            address=self.address,
            rpc_port=self.rpc_port,
            stream_port=self.stream_port,
        )
        return self.connection

    @property
    def space_center(self):
        if self.connection is None:
            self.connect()
        return self.connection.space_center

    def active_vessel_summary(self):
        vessel = self.space_center.active_vessel
        return {
            "ut": float(self.space_center.ut),
            "vessel": vessel.name,
            "body": vessel.orbit.body.name,
            "semi_major_axis": float(vessel.orbit.semi_major_axis),
            "eccentricity": float(vessel.orbit.eccentricity),
            "inclination": float(vessel.orbit.inclination),
        }

    def clear_maneuver_nodes(self, vessel=None):
        vessel = vessel or self.space_center.active_vessel
        for node in list(vessel.control.nodes):
            node.remove()

    def inject_plan(self, plan, vessel=None, clear_existing=False):
        vessel = vessel or self.space_center.active_vessel
        if clear_existing:
            self.clear_maneuver_nodes(vessel)
        created = []
        for node in plan.nodes:
            created_node = vessel.control.add_node(
                node.ut_seconds,
                prograde=node.prograde,
                normal=node.normal,
                radial=node.radial,
            )
            created.append(created_node)
        return created

    def dry_run_plan(self, plan):
        return plan.as_dict()
