# -*- coding: utf-8 -*-
"""
Created on Wed May 29 08:02:50 2019

@author: v.fave
"""

"""Trajectory decoding and text display helpers for pykep mga_1dsm genes."""

from numpy import dot,cross,linalg,around,rad2deg, array
import logging
from pykep import epoch
from scipy.linalg import norm
from pykep.core import DAY2SEC, lambert_problem, propagate_lagrangian, ic2par
from pykep.core import fb_vout
from math import sqrt, pi, cos, sin, acos, log, atan, atan2, asin, copysign
from math import acosh, sinh
from planet_packs import PACKS
import numpy as np


logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


def _safe_acos(value):
    """Return acos with round-off protection for normalized geometry."""
    return acos(float(np.clip(value, -1.0, 1.0)))


def Kdate(time,planet_pack):
    if planet_pack == "Vanilla":
        return("Y"+str(int(time//426)+1)+" D"+str(int(time%426)+1))
    elif planet_pack == "JNSQ":
        return("Y"+str(int(time//365)+1)+" D"+str(int(time%365)+1))
        
def Ktime(time,planet_pack):
    if planet_pack == "Vanilla":    
        if int(round(time//426,0)) < 1:
            return(str(int(round(time%426,0)))+"D")
        else :
            return(str(int(round(time//426,0)))+"Y "+str(int(round(time%426,0)))+"D")
    elif planet_pack == "JNSQ":
        if int(round(time//365,0)) < 1:
            return(str(int(round(time%365,0)))+"D")
        else :
            return(str(int(round(time//365,0)))+"Y "+str(int(round(time%365,0)))+"D")  


def _fmt(value, digits=3):
    return str(round(float(value), digits))
    
 
    
def _body_data(planet_pack):
    if planet_pack not in PACKS:
        raise ValueError("Unknown planet pack: " + str(planet_pack))
    pack = PACKS[planet_pack]
    return pack.body, pack.row, pack.col


def _starting_soi_radius(planet, planet_pack):
    body, row, col = _body_data(planet_pack)
    return body[row[planet.name], col['R_soi (km)']] * 1000


def _seconds_to_kdays(seconds, planet_pack):
    return float(seconds) / DAY2SEC * PACKS[planet_pack].EDY_TO_KDY


def _hyperbolic_periapsis_to_radius_time(body, vinf_speed, periapsis_radius,
                                         boundary_radius):
    """Return time from hyperbolic periapsis to a radius, in seconds."""
    vinf_speed = float(vinf_speed)
    periapsis_radius = float(periapsis_radius)
    boundary_radius = float(boundary_radius)
    if vinf_speed <= 1e-15 or periapsis_radius <= 0.0:
        return None
    if boundary_radius <= periapsis_radius:
        return 0.0
    semi_major_abs = float(body.mu_self) / (vinf_speed * vinf_speed)
    eccentricity = 1.0 + periapsis_radius / semi_major_abs
    cosh_anomaly = (boundary_radius / semi_major_abs + 1.0) / eccentricity
    if cosh_anomaly < 1.0:
        if cosh_anomaly > 1.0 - 1e-12:
            cosh_anomaly = 1.0
        else:
            return None
    anomaly = acosh(cosh_anomaly)
    return sqrt(semi_major_abs ** 3 / float(body.mu_self)) * (
        eccentricity * sinh(anomaly) - anomaly
    )


def _hyperbolic_periapsis_to_boundary_time(body, boundary_speed,
                                           periapsis_radius, boundary_radius):
    """Return time to a boundary when speed at that boundary is known."""
    boundary_speed = float(boundary_speed)
    vinf_squared = (
        boundary_speed * boundary_speed
        - 2.0 * float(body.mu_self) / float(boundary_radius)
    )
    if vinf_squared <= 0.0:
        return None
    return _hyperbolic_periapsis_to_radius_time(
        body, sqrt(vinf_squared), periapsis_radius, boundary_radius
    )


def _patched_conic_timing_estimate(first_body, target_body, departure_vinf,
                                   arrival_vinf, departure_altitude,
                                   target_periapsis_radius, planet_pack):
    """Estimate KSP-visible SOI/Pe timing offsets around a heliocentric leg."""
    departure_time = _hyperbolic_periapsis_to_boundary_time(
        first_body,
        departure_vinf,
        float(first_body.radius) + float(departure_altitude),
        _starting_soi_radius(first_body, planet_pack),
    )
    arrival_time = _hyperbolic_periapsis_to_radius_time(
        target_body,
        arrival_vinf,
        target_periapsis_radius,
        _starting_soi_radius(target_body, planet_pack),
    )
    departure_kdays = (
        _seconds_to_kdays(departure_time, planet_pack)
        if departure_time is not None else None
    )
    arrival_kdays = (
        _seconds_to_kdays(arrival_time, planet_pack)
        if arrival_time is not None else None
    )
    return departure_kdays, arrival_kdays


def _max_revs(self):
    return getattr(self, "max_revs", 0)


def _propagate_lagrangian(r, v, tof, mu):
    return propagate_lagrangian([r, v], tof, mu)


def _lambert_problem(r0, r1, tof, mu, cw=False, max_revs=0):
    return lambert_problem(r0, r1, tof, mu, cw=cw, multi_revs=max_revs)


def _lambert_v1(lambert):
    return lambert.v0[0]


def _lambert_v2(lambert):
    return lambert.v1[0]


def _ic2par(r, v, mu):
    return ic2par([r, v], mu)


     

def EJ_BurnDV(self,v_HEV,ParkingAlt,Rsoi,ejectionInclination,Alt):        
        
    Planet = self._seq[0]            
    '''
    Okay so in short, we can get the cos(Delta i)'s from the v_HEV's
    Then do the calculation in a proper way. So... cos = adj/hyp
    however, that is not the whole story. The Delta i is not the angle at SOI,
    but the angle between the hyperbolic plane the initial parking orbit plane,
    knowing that the planes intersect at Pe.
    '''
        
    v1 = np.linalg.norm(v_HEV)
    v0 = sqrt(Planet.mu_self / (Planet.radius+Alt) )
    v1 = np.sqrt(v1**2 + 2 * Planet.mu_self / (Planet.radius+Alt) -2 * Planet.mu_self / Rsoi)      
    BurnDV = np.sqrt(v0**2 + v1**2 - 2*v1*v0*cos(ejectionInclination))
    BurnDV_xy = v1-v0 
    #just to ward of occasional rounding errors leading to negative value in the root
    #will still crash in case of a real anomalous result
    if abs(BurnDV-BurnDV_xy) < 1e-8 :
        BurnDV_z = 0
    else :
        BurnDV_z = copysign(np.sqrt(BurnDV**2-BurnDV_xy**2),v_HEV[2])
    
    return BurnDV,BurnDV_xy,BurnDV_z

def EJ_Pe_Direction(self,VSOI,theta):
    if abs(VSOI[1]) < 1e-12:
        logger.warning("VSOI[1] is too small, set to 1e-12")
        VSOI[1] = 1e-12
    cosTheta = cos(theta)
    g = -VSOI[0]/VSOI[1]
    a = 1+g*g
    b = 2*g*cosTheta/VSOI[1]
    c = cosTheta**2 / (VSOI[1]**2) -1 
    if b < 0 :
        q = -0.5*(b - sqrt(b**2 - 4*a*c))
    else :
        q = -0.5*(b + sqrt(b**2 - 4*a*c))
    vx = q/a
    vy = g*vx+ cosTheta / VSOI[1]
    if cross([vx,vy,0],VSOI)[2]<0:
        vx = c/q
        vy = g*vx + cosTheta/VSOI[1]
    return np.array([vx,vy,0])
    
def EJ_angle_from_Pe(self,EJDV_at_SOI,Rsoi,Alt):
    Planet = self._seq[0]
    orbital_speed = sqrt(Planet.mu_self/(Planet.radius+Alt))
    v1 = sqrt(EJDV_at_SOI*EJDV_at_SOI +2*orbital_speed*orbital_speed-2*Planet.mu_self/Rsoi)
    e = (Planet.radius+Alt)*v1*v1/Planet.mu_self-1
    a = (Planet.radius+Alt)/(1-e)
    theta = _safe_acos((a * (1 - e * e) - Rsoi) / (e * Rsoi))
    return theta + asin(v1 * (Planet.radius+Alt) / (EJDV_at_SOI * Rsoi))

def EJ_angle_to_Prograde(self,Pe,Pro):

    pro_zeroZ = np.array([Pro[0],Pro[1],0])
    pro_zeroZ = pro_zeroZ/linalg.norm(pro_zeroZ)
    if cross(Pe,pro_zeroZ)[2] < 0 :
        return 2*pi - _safe_acos(dot(Pe,pro_zeroZ))
    else :
        return _safe_acos(dot(Pe,pro_zeroZ))

def EJ_details(self,v_EJDV_at_SOI,Rsoi,prograde_uv,Alt):
    EJDV_at_SOI = linalg.norm(v_EJDV_at_SOI)
    theta = EJ_angle_from_Pe(self,EJDV_at_SOI,Rsoi,Alt)
    ejectionDirection = v_EJDV_at_SOI/ EJDV_at_SOI

    if (abs(sin(theta)) < abs(ejectionDirection[2])) :
        logger.warning(
            "Invalid ejection geometry: theta=%s and ejdr[2]=%s",
            round(theta, 1),
            round(ejectionDirection[2], 1),
        )
        ejectionDeltaV = 99999,99999,99999
        return [ejectionDeltaV, float(90), float(90)]
    else :
        periapsisDirection = EJ_Pe_Direction(self,ejectionDirection, theta)
        ejectionAngle = EJ_angle_to_Prograde(self,periapsisDirection, prograde_uv)
        Pe_X_Ej = cross(periapsisDirection, ejectionDirection)
        ejectionInclination = _safe_acos((Pe_X_Ej/linalg.norm(Pe_X_Ej))[2])
        ejectionInclination = copysign(ejectionInclination,pi - theta)
        ejectionInclination = copysign(ejectionInclination,ejectionDirection[2])
        ejectionDeltaV = EJ_BurnDV(self,v_EJDV_at_SOI,Alt,Rsoi,ejectionInclination,Alt)
        return ejectionDeltaV, ejectionInclination, ejectionAngle;


def _vinf_vector_from_gene(gene):
    theta = 2 * pi * gene[1]
    phi = acos(2 * gene[2] - 1) - pi / 2
    return np.asarray([
        gene[3] * cos(phi) * cos(theta),
        gene[3] * cos(phi) * sin(theta),
        gene[3] * sin(phi),
    ], dtype=float)


def _equatorial_ejection_branch(direction, cosine_asymptote, periapsis_speed,
                                 circular_speed):
    """Return the cheapest equatorial-periapsis branch for one asymptote."""
    horizontal = float(np.hypot(direction[0], direction[1]))
    if horizontal <= 1e-15 or abs(cosine_asymptote) > horizontal + 1e-12:
        return None
    longitude = atan2(direction[1], direction[0])
    offset = acos(float(np.clip(cosine_asymptote / horizontal, -1.0, 1.0)))
    candidates = []
    for periapsis_longitude in (longitude + offset, longitude - offset):
        radius_hat = np.asarray([
            cos(periapsis_longitude), sin(periapsis_longitude), 0.0
        ])
        velocity_hat = direction - dot(direction, radius_hat) * radius_hat
        velocity_norm = linalg.norm(velocity_hat)
        if velocity_norm <= 1e-15:
            continue
        velocity_hat /= velocity_norm
        circular_hat = np.asarray([-radius_hat[1], radius_hat[0], 0.0])
        burn_vector = periapsis_speed * velocity_hat - circular_speed * circular_hat
        burn = linalg.norm(burn_vector)
        plane_change = _safe_acos(dot(velocity_hat, circular_hat))
        candidates.append({
            "dv": float(burn),
            "inclination": float(plane_change),
            "periapsis_direction": radius_hat,
            "prograde_dv": float(dot(burn_vector, circular_hat)),
            "radial_dv": float(dot(burn_vector, radius_hat)),
            # KSP maneuver-node normal for the equatorial parking burn follows
            # the local burn component produced by this branch geometry. The
            # PyKEP gene latitude convention is handled when decoding vinf.
            "normal_dv": float(burn_vector[2]),
        })
    return min(candidates, key=lambda candidate: candidate["dv"]) if candidates else None


def _soi_radius_from_body(first_body):
    if hasattr(first_body, "orbital_elements"):
        semi_major_axis = float(first_body.orbital_elements[0])
    else:
        semi_major_axis = float(first_body.elements(0)[0])
    return semi_major_axis * (
        float(first_body.mu_self) / float(first_body.mu_central_body)
    ) ** 0.4


def _finite_soi_ejection_angle(first_body, soi_speed, periapsis_speed,
                               radius, rsoi):
    if soi_speed <= 1e-15:
        return None
    eccentricity = radius * periapsis_speed ** 2 / first_body.mu_self - 1.0
    if eccentricity <= 1.0:
        return None
    semi_major_axis = radius / (1.0 - eccentricity)
    true_anomaly_at_soi = _safe_acos(
        (semi_major_axis * (1.0 - eccentricity * eccentricity) - rsoi)
        / (eccentricity * rsoi)
    )
    zenith_angle_at_soi = asin(
        np.clip(periapsis_speed * radius / (soi_speed * rsoi), -1.0, 1.0)
    )
    return true_anomaly_at_soi + zenith_angle_at_soi


def _finite_soi_ejection_branch(first_body, direction, soi_speed, radius,
                                rsoi, circular_speed):
    periapsis_speed_squared = (
        soi_speed ** 2
        + 2.0 * first_body.mu_self / radius
        - 2.0 * first_body.mu_self / rsoi
    )
    if periapsis_speed_squared <= 0.0:
        return None
    periapsis_speed = sqrt(periapsis_speed_squared)
    ejection_angle = _finite_soi_ejection_angle(
        first_body, soi_speed, periapsis_speed, radius, rsoi
    )
    if ejection_angle is None:
        return None
    return _equatorial_ejection_branch(
        direction, cos(ejection_angle), periapsis_speed, circular_speed
    )


def fast_ejection_from_vinf(first_body, vinf, altitude, rsoi=None):
    """Evaluate direct and SOI-split ejections from an inertial vinf vector."""
    vinf = np.asarray(vinf, dtype=float)
    if vinf.shape != (3,):
        raise ValueError("vinf must contain three Cartesian components")
    vinf_norm = float(linalg.norm(vinf))
    if rsoi is None:
        rsoi = _soi_radius_from_body(first_body)
    rsoi = float(rsoi)
    radius = float(first_body.radius + altitude)
    circular_speed = sqrt(first_body.mu_self / radius)
    if vinf_norm <= 1e-15:
        direct = {
            "dv": float((sqrt(2.0) - 1.0) * circular_speed),
            "inclination": 0.0,
            "periapsis_direction": np.asarray([1.0, 0.0, 0.0]),
            "prograde_dv": float((sqrt(2.0) - 1.0) * circular_speed),
            "radial_dv": 0.0,
            "normal_dv": 0.0,
        }
    else:
        direct = _finite_soi_ejection_branch(
            first_body, vinf / vinf_norm, vinf_norm, radius, rsoi,
            circular_speed
        )
    if direct is None:
        direct = {
            "dv": 99999.0, "inclination": pi / 2,
            "periapsis_direction": np.asarray([1.0, 0.0, 0.0]),
            "prograde_dv": 99999.0,
            "radial_dv": 0.0,
            "normal_dv": 0.0,
        }

    planar_vinf = float(np.hypot(vinf[0], vinf[1]))
    if planar_vinf <= 1e-15:
        planar_burn = (sqrt(2.0) - 1.0) * circular_speed
        planar_prograde_dv = planar_burn
        planar_radial_dv = 0.0
        planar_normal_dv = 0.0
        planar_periapsis = np.asarray([1.0, 0.0, 0.0])
    else:
        planar_direction = np.asarray([vinf[0], vinf[1], 0.0]) / planar_vinf
        planar = _finite_soi_ejection_branch(
            first_body, planar_direction, planar_vinf, radius, rsoi,
            circular_speed
        )
        if planar is None:
            planar = {
                "dv": 99999.0,
                "prograde_dv": 99999.0,
                "radial_dv": 0.0,
                "normal_dv": 0.0,
                "periapsis_direction": np.asarray([1.0, 0.0, 0.0]),
            }
        planar_burn = planar["dv"]
        planar_prograde_dv = planar["prograde_dv"]
        planar_radial_dv = planar["radial_dv"]
        planar_normal_dv = planar["normal_dv"]
        planar_periapsis = planar["periapsis_direction"]
    soi_correction = abs(float(vinf[2]))
    split_normal_dv = float(vinf[2])
    split_dv = float(planar_burn + soi_correction)

    if direct["dv"] <= split_dv:
        strategy = "direct"
        selected_dv = direct["dv"]
        selected_inclination = direct["inclination"]
        periapsis_direction = direct["periapsis_direction"]
        selected_normal_dv = direct["normal_dv"]
        selected_prograde_dv = direct["prograde_dv"]
        selected_radial_dv = direct["radial_dv"]
        parking_node_dv = direct["dv"]
        parking_node_prograde_dv = direct["prograde_dv"]
        parking_node_radial_dv = direct["radial_dv"]
        parking_node_normal_dv = direct["normal_dv"]
        soi_normal_correction_dv = 0.0
        selected_normal_dv_location = "parking_orbit"
    else:
        strategy = "split_soi"
        selected_dv = split_dv
        selected_inclination = 0.0
        periapsis_direction = planar_periapsis
        selected_normal_dv = split_normal_dv
        selected_prograde_dv = planar_prograde_dv
        selected_radial_dv = planar_radial_dv
        parking_node_dv = planar_burn
        parking_node_prograde_dv = planar_prograde_dv
        parking_node_radial_dv = planar_radial_dv
        parking_node_normal_dv = planar_normal_dv
        soi_normal_correction_dv = split_normal_dv
        selected_normal_dv_location = "soi"
    return {
        "dv": float(selected_dv),
        "strategy": strategy,
        "direct_dv": float(direct["dv"]),
        "direct_prograde_dv": float(direct["prograde_dv"]),
        "direct_radial_dv": float(direct["radial_dv"]),
        "direct_normal_dv": float(direct["normal_dv"]),
        "direct_inclination": float(direct["inclination"]),
        "direct_periapsis_direction": direct["periapsis_direction"],
        "planar_dv": float(planar_burn),
        "planar_prograde_dv": float(planar_prograde_dv),
        "planar_radial_dv": float(planar_radial_dv),
        "planar_normal_dv": float(planar_normal_dv),
        "planar_periapsis_direction": planar_periapsis,
        "soi_correction_dv": float(soi_correction),
        "split_normal_dv": float(split_normal_dv),
        "split_dv": split_dv,
        "selected_prograde_dv": float(selected_prograde_dv),
        "selected_radial_dv": float(selected_radial_dv),
        "selected_normal_dv": float(selected_normal_dv),
        "parking_node_dv": float(parking_node_dv),
        "parking_node_prograde_dv": float(parking_node_prograde_dv),
        "parking_node_radial_dv": float(parking_node_radial_dv),
        "parking_node_normal_dv": float(parking_node_normal_dv),
        "soi_normal_correction_dv": float(soi_normal_correction_dv),
        "selected_normal_dv_location": selected_normal_dv_location,
        "inclination": float(selected_inclination),
        "periapsis_direction": periapsis_direction,
        "vinf": vinf,
    }


def fast_ejection_from_gene(first_body, gene, altitude, rsoi=None):
    """Evaluate canonical ejection after decoding a PyKEP MGA gene."""
    return fast_ejection_from_vinf(
        first_body,
        _vinf_vector_from_gene(gene),
        altitude,
        rsoi,
    )


def ejection_from_gene(self, gene, rsoi, altitude):
    """Return the canonical parking-orbit ejection calculation for a gene."""
    if not hasattr(self, "_decode_times_and_vinf"):
        from pykep.trajopt import mga_1dsm
        extracted = self.extract(mga_1dsm)
        if extracted is None:
            raise TypeError("Expected an mga_1dsm problem for ejection decoding")
        self = extracted
    result = fast_ejection_from_gene(self._seq[0], gene, altitude, rsoi)
    _, departure_velocity = self._seq[0].eph(epoch(gene[0]))
    prograde = departure_velocity / linalg.norm(departure_velocity)
    result["angle"] = float(EJ_angle_to_Prograde(
        self, result["periapsis_direction"], prograde
    ))
    result["full_dv"] = result["dv"]
    result["vertical_dv"] = result["selected_normal_dv"]
    return result
def transx(self,x,alt=100000,planet_pack="Vanilla"):
    if planet_pack == "Vanilla":
        Edy2Kdy = 4
    elif planet_pack == "JNSQ":
        Edy2Kdy = 2
    # 1 -  we 'decode' the chromosome recording the various times of flight
    # (days) in the list T and the cartesian components of vinf
    # n_legs = len(self._seq)

    T, Vinfx, Vinfy, Vinfz = self._decode_times_and_vinf(x)
    t_P = list([None] * (self.n_legs+1))
    r_P = list([None] * (self.n_legs+1))
    v_P = list([None] * (self.n_legs+1))
    DV  = list([None] * (self.n_legs+1))
        
    for i, planet in enumerate(self._seq):
        t_P[i] = epoch(x[0] + sum(T[0:i]))
        # Here one issue is that kerbin epoch zero is Y1D1. We have to add 1Kdy to the print result (BUT NOT TO THE CALC)
        r_P[i], v_P[i] = self._seq[i].eph(t_P[i]) 
        
    fward_P = list([None] * (self.n_legs+1))
    plane_P = list([None] * (self.n_legs+1))
    oward_P = list([None] * (self.n_legs+1))
    Vinf = [Vinfx,Vinfy,Vinfz]

    for i, planet in enumerate(self._seq):
        fward_P[i] = v_P[i] / linalg.norm(v_P[i])
        plane_P[i] = cross(v_P[i], r_P[i])
        plane_P[i] = plane_P[i] / linalg.norm(plane_P[i])
        oward_P[i] = cross(plane_P[i], fward_P[i])

    # 3 - We start with the first leg
    Rsoi = _starting_soi_radius(self._seq[0], planet_pack)

    ejection = ejection_from_gene(self, x, Rsoi, alt)
    DV_Ej = (
        ejection["full_dv"], ejection["planar_dv"], ejection["vertical_dv"]
    )
    DV_i = ejection["inclination"] * 180/pi
    ejection_angle = ejection["angle"]
    DV_Sum = ejection["dv"]
        
    pe = alt+self._seq[0].radius
    a = -1*self._seq[0].mu_self/(x[3]*x[3])
    e = (a-pe)/a

    Eta = _safe_acos(-1/e)
    Eta_Rsoi = _safe_acos((a * (1 - e * e) - Rsoi) / (e * Rsoi));
    Delta = asin(1/e)*2
    Gamma = atan(dot(r_P[0],Vinf)/(norm(Vinf)*norm(r_P[0])))
    
    v0 = [a + b for a, b in zip(v_P[0], [Vinfx, Vinfy, Vinfz])]
    
    k_el = _ic2par(r_P[0],v0,self.common_mu) 
    Ap = k_el[0]*(1+k_el[1])
    Pe = k_el[0]*(1-k_el[1])
       
    print("")
    print("First Leg:               " + self._seq[0].name + " to " + self._seq[1].name)
    print("Departure:               " + str(round(t_P[0].mjd2000*Edy2Kdy,1)) + " KUT ("+Kdate(round(t_P[0].mjd2000*Edy2Kdy,1),planet_pack)+")")
    print("Departure precise:       " + _fmt(t_P[0].mjd2000*Edy2Kdy,4) + " KUT")
    print("Duration:                " + Ktime(T[0]*Edy2Kdy,planet_pack))
    print("Duration precise:        " + _fmt(T[0]*Edy2Kdy,4) + " days")
    print("VINF:                    " + str(round(x[3],1)) + " m/s")
    print("VINF (x,y,z):            " + str(around(Vinf,1))+ " m/s")
    #Rp and Vp are quite useless. Printing the pe/ap would make a lot more sense. Ap = SMA*(1+e) & Pe = SMA*(1-e)
    print("Ap:                      " + str(int(round(Ap/1000,0)))+ " km")
    print("Pe:                      " + str(int(round(Pe/1000,0)))+ " km")
    print("Ejection DV:             " + str(round(DV_Ej[0],1)) + " m/s, from "+str(round(alt/1000,1))+ " km parking orbit")     
    z_location = "SOI" if ejection["selected_normal_dv_location"] == "soi" else "parking orbit"
    print("Ejection DV xy and z:    " + str(round(DV_Ej[1],1)) + " m/s (xy) and "+ str(round(DV_Ej[2],1)) + " m/s (selected z at " + z_location + ")")
    print("Parking node components: " + str(round(ejection["parking_node_prograde_dv"],1)) + " m/s prograde, " + str(round(ejection["parking_node_normal_dv"],1)) + " m/s normal, " + str(round(ejection["parking_node_radial_dv"],1)) + " m/s radial")
    print("Parking node precise:    " + _fmt(ejection["parking_node_prograde_dv"]) + " m/s prograde, " + _fmt(ejection["parking_node_normal_dv"]) + " m/s normal, " + _fmt(ejection["parking_node_radial_dv"]) + " m/s radial")
    print("Ejection angle precise:  " + _fmt(ejection_angle*180/pi,3) + " degrees from prograde")
    print("Ejection angle:          " + str(round(ejection_angle*180/pi,1)) + " ° from prograde")
    print("Ejection inclination:    " + str(round(DV_i,1)) + " °")   
      
    planar_vinf = array([Vinfx, Vinfy, 0.0])
    soi_correction_vector = array([0.0, 0.0, Vinfz])
    soi_correction_time = _hyperbolic_periapsis_to_boundary_time(
        self._seq[0],
        norm(planar_vinf),
        alt + self._seq[0].radius,
        Rsoi,
    )
    soi_correction_components = None
    if soi_correction_time is not None:
        soi_r, soi_v = _propagate_lagrangian(
            r_P[0],
            array(v_P[0]) + planar_vinf,
            soi_correction_time,
            self.common_mu,
        )
        soi_prograde_uv = array(soi_v) / linalg.norm(soi_v)
        soi_plane_uv = cross(soi_v, soi_r)
        soi_plane_uv = soi_plane_uv / linalg.norm(soi_plane_uv)
        soi_radial_uv = cross(soi_plane_uv, soi_prograde_uv)
        soi_correction_components = (
            dot(soi_prograde_uv, soi_correction_vector),
            -1 * dot(soi_plane_uv, soi_correction_vector),
            dot(soi_radial_uv, soi_correction_vector),
            _seconds_to_kdays(soi_correction_time, planet_pack),
        )

    if ejection["strategy"] == "split_soi":
        print("Split burn is optimal :  " + str(round(DV_Ej[1],1)) + " m/s (xy) at Pe and "+str(round(ejection["split_normal_dv"],1))+ " m/s (z) at SOI")
        print("SOI normal correction:   " + str(round(ejection["soi_normal_correction_dv"],1)) + " m/s")
        if soi_correction_components is not None:
            print("SOI correction precise:  T+" + _fmt(soi_correction_components[3],4) + " days, " + _fmt(t_P[0].mjd2000*Edy2Kdy + soi_correction_components[3],4) + " KUT")
            print("  node components:       " + _fmt(soi_correction_components[0]) + " m/s prograde, " + _fmt(soi_correction_components[1]) + " m/s normal, " + _fmt(soi_correction_components[2]) + " m/s radial")
    else:
        print("Direct inclined burn is optimal: " + str(round(ejection["direct_normal_dv"],1)) + " m/s (z) from parking orbit")
        print("Split-SOI alternative:   " + _fmt(ejection["planar_prograde_dv"]) + " m/s prograde, " + _fmt(ejection["planar_normal_dv"]) + " m/s normal, " + _fmt(ejection["planar_radial_dv"]) + " m/s radial at parking")
        print("  plus SOI normal:       " + _fmt(ejection["split_normal_dv"]) + " m/s (total +" + _fmt(ejection["split_dv"] - ejection["direct_dv"]) + " m/s vs direct)")
        if soi_correction_components is not None:
            print("  SOI node estimate:     T+" + _fmt(soi_correction_components[3],4) + " days, " + _fmt(t_P[0].mjd2000*Edy2Kdy + soi_correction_components[3],4) + " KUT")
            print("    components:          " + _fmt(soi_correction_components[0]) + " m/s prograde, " + _fmt(soi_correction_components[1]) + " m/s normal, " + _fmt(soi_correction_components[2]) + " m/s radial")

    r, v = _propagate_lagrangian(
        r_P[0], v0, x[4] * T[0] * DAY2SEC, self.common_mu)

    first_dsm_met = x[4] * T[0] * Edy2Kdy
    first_dsm_epoch = t_P[0].mjd2000 * Edy2Kdy + first_dsm_met
    print("DSM after :              " + Ktime(first_dsm_met,planet_pack))
    print("DSM epoch precise:       T+" + _fmt(first_dsm_met,4) + " days, " + _fmt(first_dsm_epoch,4) + " KUT")

    # Lambert arc to reach seq[1]
    dt = (1 - x[4]) * T[0] * DAY2SEC
    l = _lambert_problem(
        r, r_P[1], dt, self.common_mu, cw = False, max_revs=_max_revs(self))
    v_end_l = _lambert_v2(l)
    v_beg_l = _lambert_v1(l)

    # First DSM occuring at time nu1*T1
    DV[0] = norm([a - b for a, b in zip(v_beg_l, v)])
    DV_vect = [a - b for a, b in zip(v_beg_l, v)]
    # repère local
    progrd_uv   = array(v) / linalg.norm(v)
    plane_uv   = cross(v, r)
    plane_uv   = plane_uv / linalg.norm(plane_uv)
    radial_uv   = cross(plane_uv, progrd_uv)
    print("DSM magnitude:           " + str(round(DV[0],1)) + " m/s")
    print("Prograde:                " + str(np.round(dot(progrd_uv, DV_vect), 1)) + " m/s")
    print("Normal:                  " + str(-1*np.round(dot(plane_uv , DV_vect), 1)) + " m/s") 
    print("Radial:                  " + str(np.round(dot(radial_uv, DV_vect), 1)) + " m/s")
    print("DSM precise:             " + _fmt(dot(progrd_uv, DV_vect)) + " m/s prograde, " + _fmt(-1*dot(plane_uv, DV_vect)) + " m/s normal, " + _fmt(dot(radial_uv, DV_vect)) + " m/s radial")
    try:
        from _KSPTranslator import format_translation_summary
        from _KSPTranslator import translate_first_leg
        translation = translate_first_leg(
            self, x, planet_pack=planet_pack, parking_altitude=alt
        )
        for line in format_translation_summary(translation, planet_pack=planet_pack):
            print(line)
    except Exception as exc:
        logger.warning("KSP translation diagnostics failed: %s", exc)
    DV_Sum += DV[0]
    
    # 4 - And we proceed with each successive leg
    for i in range(1, self.n_legs):
        print("\nleg no. " + str(i + 1) + ":               " +
                  self._seq[i].name + " to " + self._seq[i + 1].name)
        print("Duration:                " + Ktime(T[i]*Edy2Kdy,planet_pack))
        # Fly-by
        v_out = fb_vout(v_end_l, v_P[i], x[
                7 + (i - 1) * 4] * self._seq[i].radius, x[6 + (i - 1) * 4], self._seq[i].mu_self)
        incoming_vinf_vector = array(v_end_l) - array(v_P[i])
        outgoing_vinf_vector = array(v_out) - array(v_P[i])
        incoming_hat = incoming_vinf_vector / linalg.norm(incoming_vinf_vector)
        outgoing_hat = outgoing_vinf_vector / linalg.norm(outgoing_vinf_vector)
        flyby_plane_normal = cross(incoming_hat, outgoing_hat)
        flyby_plane_normal = flyby_plane_normal / linalg.norm(flyby_plane_normal)
        hyperbola_inclination = acos(
            float(np.clip(flyby_plane_normal[2], -1.0, 1.0))
        ) * 180 / pi
        print(
            "Fly-by epoch:            " + Kdate(t_P[i].mjd2000*Edy2Kdy,planet_pack)+" ("+Ktime(t_P[i].mjd2000*Edy2Kdy-t_P[0].mjd2000*Edy2Kdy,planet_pack)+" after launch)")
        print(
            "Fly-by radius:           " + str(round(x[7 + (i - 1) * 4],2)) + " planetary radii ("+str(round((x[7 + (i - 1) * 4]-1)*self._seq[i].radius/1000.0,1))+" km altitude)")
        print(
            "Beta plane angle:        " + str(round(x[6 + (i - 1) * 4]*180/pi,0))+"° ")
        print("Hyperbola inclination:   " + str(round(hyperbola_inclination, 2))+" deg")
        print(
            "Flyby plane normal:      ["
            + ", ".join(str(round(v, 4)) for v in flyby_plane_normal)
            + "]"
        )
        # s/c propagation before the DSM
        r, v = _propagate_lagrangian(
            r_P[i], v_out, x[8 + (i - 1) * 4] * T[i] * DAY2SEC, self.common_mu)
        print("DSM after :              " + Ktime(x[8 + (i - 1) * 4] * T[i]*Edy2Kdy,planet_pack) + " after flyby ("+Ktime((x[8 + (i - 1) * 4] * T[i]+t_P[i].mjd2000-t_P[0].mjd2000)*Edy2Kdy,planet_pack)+" days after launch)")
        # Lambert arc to reach Earth during (1-nu2)*T2 (second segment)
        dt = (1 - x[8 + (i - 1) * 4]) * T[i] * DAY2SEC
        l = _lambert_problem(r, r_P[i + 1], dt, self.common_mu, cw=False, max_revs=_max_revs(self))
        v_end_l = _lambert_v2(l)
        v_beg_l = _lambert_v1(l)
        # DSM occuring at time nu2*T2
        # DV in vector form, decompose the node components
        DV_vect = [a - b for a, b in zip(v_beg_l, v)]
        # repère local
        progrd_uv   = array(v) / linalg.norm(v)
        plane_uv   = cross(v, r)
        plane_uv   = plane_uv / linalg.norm(plane_uv)
        radial_uv   = cross(plane_uv, progrd_uv)
        
        DV[i] = norm([a - b for a, b in zip(v_beg_l, v)])
        
        DV_Sum += DV[i]
        print("DSM magnitude:           " + str(round(DV[i],1)) + " m/s")
        print("Prograde:                " + str(np.round(dot(progrd_uv, DV_vect), 1)) + " m/s")
        print("Normal:                  " + str(-1*np.round(dot(plane_uv , DV_vect), 1)) + " m/s")
        print("Radial:                  " + str(np.round(dot(radial_uv, DV_vect), 1)) + " m/s")


    DV[-1] = norm([a - b for a, b in zip(v_end_l, v_P[-1])])

    if self._orbit_insertion:
                # In this case we compute the insertion DV as a single pericenter
                # burn
        print(DV[-1])
        DVper = np.sqrt(DV[-1] * DV[-1] + 2 *
                                self._seq[-1].mu_self / self._rp_target)
        DVper2 = np.sqrt(2 * self._seq[-1].mu_self / self._rp_target -
                                self._seq[-1].mu_self / self._rp_target * (1. - self._e_target))
        #DVper = np.sqrt(DV[-1] * DV[-1] + 2 * self._seq[-1].mu_self / (self._rp_target) - 2 * self._seq[-1].mu_self / Rsoi_tgt)
        #DVper2 = np.sqrt(self._seq[-1].mu_self / (self._rp_target) * (1. - self._e_target))
        DV_Inj = np.abs(DVper - DVper2)   
        

        #v0 = sqrt(self._seq[-1].mu_self / (self._seq[-1].radius+350000) )
        #v1 = np.sqrt(DV[-1]**2 + 2 * self._seq[-1].mu_self / (self._seq[-1].radius+350000) -2 * self._seq[-1].mu_self / Rsoi) - v0
        

    arrival_met = t_P[-1].mjd2000*Edy2Kdy - t_P[0].mjd2000*Edy2Kdy
    target_periapsis_radius = (
        self._rp_target if self._orbit_insertion else self._seq[-1].radius
    )
    departure_soi_kdays, arrival_soi_to_pe_kdays = (
        _patched_conic_timing_estimate(
            self._seq[0],
            self._seq[-1],
            norm(Vinf),
            DV[-1],
            alt,
            target_periapsis_radius,
            planet_pack,
        )
    )
    print("\nArrival at " + self._seq[-1].name)
    print("Arrival epoch:           " + str(round(t_P[-1].mjd2000*Edy2Kdy,1)) + " KUT ("+Kdate(t_P[-1].mjd2000*Edy2Kdy,planet_pack)+", T+"+Ktime(arrival_met,planet_pack)+" after launch, heliocentric)")
    print("Total mission time:      " + str(round(sum(T)*Edy2Kdy,1)) + " days ("+Ktime(sum(T)*Edy2Kdy,planet_pack)+")")   
    if departure_soi_kdays is not None and arrival_soi_to_pe_kdays is not None:
        print("SOI timing diagnostics:  " + str(round(departure_soi_kdays,1)) + " days departure Pe->SOI, " + str(round(arrival_soi_to_pe_kdays,1)) + " days final SOI->Pe")
        print("  KSP node TOF note:     compare the game Pe time to the heliocentric T+" + Ktime(arrival_met,planet_pack) + " above")
    print("Arrival Vinf:            " + str(round(DV[-1],1)) + " m/s")     
    print("Total DV w.o. IJB :      " + str(round(DV_Sum,1)) + " m/s")    

    if self._orbit_insertion:
        DV_Sum += DV_Inj
        print("Target pe is :           "+str(round(self._rp_target/1000,1))+" km with e = "+str(round(self._e_target,2)))    
        print("Injection DV:            " + str(round(DV_Inj,1)) + " m/s")  
        print("Total DV with IJB :      " + str(round(DV_Sum,1)) + " m/s")  
    elif self._add_vinf_arr  :
        DV_Sum += DV[-1]
        print("Total DV with V_inf_arr :  " + str(round(DV_Sum,1)) + " m/s")
        
        
def decode_trajectory(self,x,alt=100000,planet_pack="Vanilla"):
    body, row, col = _body_data(planet_pack)
    Edy2Kdy = PACKS[planet_pack].EDY_TO_KDY
        
    ''' 1 -  we 'decode' the chromosome recording the various times of flight (days) 
    in the list T and the cartesian components of vinf'''
    n_legs = len(self._seq)

    T, Vinfx, Vinfy, Vinfz = self._decode_times_and_vinf(x)
    Vinf = [Vinfx,Vinfy,Vinfz]
    
    t_P = list([None] * (n_legs))
    r_P = list([None] * (n_legs))
    v_P = list([None] * (n_legs))
    DV = list([None] * (n_legs))
        
    for i, planet in enumerate(self._seq):
        t_P[i] = epoch(x[0] + sum(T[0:i]))
        r_P[i], v_P[i] = self._seq[i].eph(t_P[i]) 
        
    fward_P = list([None] * (n_legs))

    for i, planet in enumerate(self._seq):
        fward_P[i] = v_P[i] / linalg.norm(v_P[i])
     
    ''' we compute the ejection DV accounting for inclination and possible 
    split burn if more efficient'''
    
    Rsoi = body[row[self._seq[0].name],col['R_soi (km)']]*1000
    ejection = ejection_from_gene(self, x, Rsoi, alt)
    DV_Ej = (
        ejection["full_dv"], ejection["planar_dv"], ejection["vertical_dv"]
    )
    ejection_inclination = ejection["inclination"]
    ejection_angle = ejection["angle"]
    ejection_dv = ejection["dv"]

    DV_Sum = ejection_dv

    v0 = [a + b for a, b in zip(v_P[0], [Vinfx, Vinfy, Vinfz])]
    r, v = _propagate_lagrangian(
        r_P[0], v0, x[4] * T[0] * DAY2SEC, self.common_mu)


    '''Lambert arc to reach seq[1]'''
    dt = (1 - x[4]) * T[0] * DAY2SEC
    l = _lambert_problem(r, r_P[1], dt, self.common_mu, cw = False, max_revs=_max_revs(self))
    v_end_l = _lambert_v2(l)
    v_beg_l = _lambert_v1(l)

    '''First DSM occuring at time nu1*T1'''
    DV[0] = norm([a - b for a, b in zip(v_beg_l, v)])

    DV_Sum += DV[0]
    
    '''And we proceed with each successive leg'''
    for i in range(1, self.n_legs):
        # Fly-by
        v_out = fb_vout(v_end_l, v_P[i], x[7 + (i-1) * 4] * self._seq[i].radius, x[6 + (i-1) * 4], self._seq[i].mu_self)
        # s/c propagation before the DSM
        r, v = _propagate_lagrangian(
            r_P[i], v_out, x[8 + (i-1) * 4] * T[i] * DAY2SEC, self.common_mu)
        # Lambert arc to reach Earth during (1-nu2)*T2 (second segment)
        dt = (1 - x[8 + (i-1) * 4]) * T[i] * DAY2SEC
        l = _lambert_problem(r, r_P[i + 1], dt, self.common_mu, cw=False, max_revs=_max_revs(self))
        v_end_l = _lambert_v2(l)
        v_beg_l = _lambert_v1(l)
        # DSM occuring at time nu2*T2
        DV[i] = norm([a - b for a, b in zip(v_beg_l, v)])
        DV_Sum += DV[i]   

    DV[-1] = norm([a - b for a, b in zip(v_end_l, v_P[-1])])
    DV_Inj = None
    target_periapsis_radius = self._seq[-1].radius
    if self._orbit_insertion:
        # In this case we compute the insertion DV as a single pericenter
        # burn. We assume that the final orbit is coplanar with the incoming hyperbola plane.
        target_periapsis_radius = self._rp_target
        DVper = np.sqrt(DV[-1] * DV[-1] + 2 *
                                 self._seq[-1].mu_self / self._rp_target)
        DVper2 = np.sqrt(2 * self._seq[-1].mu_self / self._rp_target -
                                 self._seq[-1].mu_self / self._rp_target * (1. - self._e_target))
        DV_Inj = np.abs(DVper - DVper2)    
            
    dsm_dv_total = sum(DV[:-1])
    dv_without_arrival = ejection_dv + dsm_dv_total
    dv_with_arrival_vinf = dv_without_arrival + DV[-1]
    dv_with_capture = None
    
    if self._orbit_insertion:
        DV_Sum += DV_Inj
        dv_with_capture = DV_Sum
    elif self._add_vinf_arr  :
        DV_Sum += DV[-1]

    departure_soi_kdays, arrival_soi_to_pe_kdays = (
        _patched_conic_timing_estimate(
            self._seq[0],
            self._seq[-1],
            norm(Vinf),
            DV[-1],
            alt,
            target_periapsis_radius,
            planet_pack,
        )
    )
    heliocentric_arrival_met = sum(T) * Edy2Kdy
    ksp_translation = None
    try:
        from _KSPTranslator import translate_first_leg
        ksp_translation = translate_first_leg(
            self, x, planet_pack=planet_pack, parking_altitude=alt
        ).as_dict()
    except Exception as exc:
        logger.warning("KSP translation diagnostics failed: %s", exc)

    return {
        "objective_dv": DV_Sum,
        "ejection_dv": ejection_dv,
        "ejection_strategy": ejection["strategy"],
        "ejection_direct_dv": ejection["direct_dv"],
        "ejection_direct_dv_prograde": ejection["direct_prograde_dv"],
        "ejection_direct_dv_radial": ejection["direct_radial_dv"],
        "ejection_direct_dv_z": ejection["direct_normal_dv"],
        "ejection_split_dv": ejection["split_dv"],
        "ejection_split_dv_z": ejection["split_normal_dv"],
        "ejection_dv_prograde": ejection["selected_prograde_dv"],
        "ejection_dv_radial": ejection["selected_radial_dv"],
        "ejection_soi_correction_dv": ejection["soi_correction_dv"],
        "ejection_dv_full": DV_Ej[0],
        "ejection_dv_xy": DV_Ej[1],
        "ejection_dv_z": DV_Ej[2],
        "ejection_dv_z_location": ejection["selected_normal_dv_location"],
        "parking_node_dv": ejection["parking_node_dv"],
        "parking_node_prograde_dv": ejection["parking_node_prograde_dv"],
        "parking_node_radial_dv": ejection["parking_node_radial_dv"],
        "parking_node_normal_dv": ejection["parking_node_normal_dv"],
        "soi_normal_correction_dv": ejection["soi_normal_correction_dv"],
        "ejection_inclination": ejection_inclination,
        "ejection_angle": ejection_angle,
        "dsm_dv": DV[:-1],
        "dsm_dv_total": dsm_dv_total,
        "arrival_vinf": DV[-1],
        "capture_dv": DV_Inj,
        "dv_without_arrival": dv_without_arrival,
        "dv_with_arrival_vinf": dv_with_arrival_vinf,
        "dv_with_capture": dv_with_capture,
        "t0": round(t_P[0].mjd2000*Edy2Kdy,1),
        "tof": round(sum(T)*Edy2Kdy,1),
        "heliocentric_arrival_met": heliocentric_arrival_met,
        "departure_soi_escape_time": departure_soi_kdays,
        "arrival_soi_to_periapsis_time": arrival_soi_to_pe_kdays,
        "arrival_periapsis_altitude_estimate": (
            target_periapsis_radius - self._seq[-1].radius
        ),
        "ksp_translation": ksp_translation,
        "ejection_vinf": norm(Vinf),
        "vinf_vector": Vinf,
        "tof_by_leg": [t*Edy2Kdy for t in T],
        "planet_pack": planet_pack,
    }


def decode_dV_tof(self,x,alt=100000,planet_pack="Vanilla"):
    decoded = decode_trajectory(self, x, alt=alt, planet_pack=planet_pack)
    return (
        decoded["objective_dv"],
        decoded["t0"],
        decoded["tof"],
        decoded["ejection_vinf"],
    )


