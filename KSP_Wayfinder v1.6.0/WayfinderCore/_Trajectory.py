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
from planet_packs import PACKS
import numpy as np


logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


def _safe_acos(value):
    """Return acos with round-off protection for normalized geometry."""
    return acos(float(np.clip(value, -1.0, 1.0)))


def Kdate(time,planet_pack):
    if planet_pack == "Vanilla":
        return("Y"+str(int(round(1+time//426,0)))+" D"+str(int(round(time%426,0))))
    elif planet_pack == "JNSQ":
        return("Y"+str(int(round(1+time//365,0)))+" D"+str(int(round(time%365,0))))
        
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
    
 
    
def _body_data(planet_pack):
    if planet_pack not in PACKS:
        raise ValueError("Unknown planet pack: " + str(planet_pack))
    pack = PACKS[planet_pack]
    return pack.body, pack.row, pack.col


def _starting_soi_radius(planet, planet_pack):
    body, row, col = _body_data(planet_pack)
    return body[row[planet.name], col['R_soi (km)']] * 1000


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
        burn = linalg.norm(
            periapsis_speed * velocity_hat - circular_speed * circular_hat
        )
        plane_change = _safe_acos(dot(velocity_hat, circular_hat))
        candidates.append({
            "dv": float(burn),
            "inclination": float(plane_change),
            "periapsis_direction": radius_hat,
        })
    return min(candidates, key=lambda candidate: candidate["dv"]) if candidates else None


def fast_ejection_from_gene(first_body, gene, altitude):
    """Evaluate direct and SOI-split ejections using closed-form geometry."""
    vinf = _vinf_vector_from_gene(gene)
    vinf_norm = float(linalg.norm(vinf))
    radius = float(first_body.radius + altitude)
    circular_speed = sqrt(first_body.mu_self / radius)
    if vinf_norm <= 1e-15:
        direct = {
            "dv": float((sqrt(2.0) - 1.0) * circular_speed),
            "inclination": 0.0,
            "periapsis_direction": np.asarray([1.0, 0.0, 0.0]),
        }
    else:
        eccentricity = 1.0 + radius * vinf_norm ** 2 / first_body.mu_self
        direct = _equatorial_ejection_branch(
            vinf / vinf_norm,
            -1.0 / eccentricity,
            sqrt(vinf_norm ** 2 + 2.0 * first_body.mu_self / radius),
            circular_speed,
        )
    if direct is None:
        direct = {
            "dv": 99999.0, "inclination": pi / 2,
            "periapsis_direction": np.asarray([1.0, 0.0, 0.0]),
        }

    planar_vinf = float(np.hypot(vinf[0], vinf[1]))
    if planar_vinf <= 1e-15:
        planar_burn = (sqrt(2.0) - 1.0) * circular_speed
        planar_periapsis = np.asarray([1.0, 0.0, 0.0])
    else:
        planar_direction = np.asarray([vinf[0], vinf[1], 0.0]) / planar_vinf
        planar_eccentricity = (
            1.0 + radius * planar_vinf ** 2 / first_body.mu_self
        )
        planar = _equatorial_ejection_branch(
            planar_direction,
            -1.0 / planar_eccentricity,
            sqrt(planar_vinf ** 2 + 2.0 * first_body.mu_self / radius),
            circular_speed,
        )
        planar_burn = planar["dv"]
        planar_periapsis = planar["periapsis_direction"]
    soi_correction = abs(float(vinf[2]))
    split_dv = float(planar_burn + soi_correction)

    if direct["dv"] <= split_dv:
        strategy = "direct"
        selected_dv = direct["dv"]
        selected_inclination = direct["inclination"]
        periapsis_direction = direct["periapsis_direction"]
    else:
        strategy = "split_soi"
        selected_dv = split_dv
        selected_inclination = 0.0
        periapsis_direction = planar_periapsis
    return {
        "dv": float(selected_dv),
        "strategy": strategy,
        "direct_dv": float(direct["dv"]),
        "direct_inclination": float(direct["inclination"]),
        "planar_dv": float(planar_burn),
        "soi_correction_dv": float(soi_correction),
        "split_dv": split_dv,
        "inclination": float(selected_inclination),
        "periapsis_direction": periapsis_direction,
        "vinf": vinf,
    }


def ejection_from_gene(self, gene, rsoi, altitude):
    """Return the canonical parking-orbit ejection calculation for a gene."""
    if not hasattr(self, "_decode_times_and_vinf"):
        from pykep.trajopt import mga_1dsm
        extracted = self.extract(mga_1dsm)
        if extracted is None:
            raise TypeError("Expected an mga_1dsm problem for ejection decoding")
        self = extracted
    result = fast_ejection_from_gene(self._seq[0], gene, altitude)
    _, departure_velocity = self._seq[0].eph(epoch(gene[0]))
    prograde = departure_velocity / linalg.norm(departure_velocity)
    result["angle"] = float(EJ_angle_to_Prograde(
        self, result["periapsis_direction"], prograde
    ))
    result["full_dv"] = result["direct_dv"]
    result["vertical_dv"] = result["soi_correction_dv"]
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

    DV_Ej,DV_i = EJ_details(self,Vinf,Rsoi,fward_P[0],alt)[0:2]
    DV_i = DV_i * 180/pi
    #4 - depending on case, split the ej burn or not.
    if DV_Ej[0] < (DV_Ej[1]+abs(Vinfz)):
        DV_Sum = DV_Ej[0]
    else :
        DV_Sum = DV_Ej[1]+abs(Vinfz)
        
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
    print("Duration:                " + Ktime(T[0]*Edy2Kdy,planet_pack))
    print("VINF:                    " + str(round(x[3],1)) + " m/s")
    print("VINF (x,y,z):            " + str(around(Vinf,1))+ " m/s")
    #Rp and Vp are quite useless. Printing the pe/ap would make a lot more sense. Ap = SMA*(1+e) & Pe = SMA*(1-e)
    print("Ap:                      " + str(int(round(Ap/1000,0)))+ " km")
    print("Pe:                      " + str(int(round(Pe/1000,0)))+ " km")
    print("Ejection DV:             " + str(round(DV_Ej[0],1)) + " m/s, from "+str(round(alt/1000,1))+ " km parking orbit")     
    print("Ejection DV xy and z:    " + str(round(DV_Ej[1],1)) + " m/s (xy) and "+ str(round(DV_Ej[2],1)) + " m/s (z)")
    print("Ejection inclination:    " + str(round(DV_i,1)) + " °")   
      
    if DV_Ej[0] > (DV_Ej[1]+abs(Vinfz)):
        print("Split burn is optimal :  " + str(round(DV_Ej[1],1)) + " m/s (xy) at Pe and "+str(round(Vinfz,1))+ " m/s (z) at SOI")   

    r, v = _propagate_lagrangian(
        r_P[0], v0, x[4] * T[0] * DAY2SEC, self.common_mu)

    print("DSM after :              " + Ktime(x[4] * T[0]*Edy2Kdy,planet_pack))

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
    DV_Sum += DV[0]
    
    # 4 - And we proceed with each successive leg
    for i in range(1, self.n_legs):
        print("\nleg no. " + str(i + 1) + ":               " +
                  self._seq[i].name + " to " + self._seq[i + 1].name)
        print("Duration:                " + Ktime(T[i]*Edy2Kdy,planet_pack))
        # Fly-by
        v_out = fb_vout(v_end_l, v_P[i], x[
                7 + (i - 1) * 4] * self._seq[i].radius, x[6 + (i - 1) * 4], self._seq[i].mu_self)
        print(
            "Fly-by epoch:            " + Kdate(t_P[i].mjd2000*Edy2Kdy,planet_pack)+" ("+Ktime(t_P[i].mjd2000*Edy2Kdy-t_P[0].mjd2000*Edy2Kdy,planet_pack)+" after launch)")
        print(
            "Fly-by radius:           " + str(round(x[7 + (i - 1) * 4],2)) + " planetary radii ("+str(round((x[7 + (i - 1) * 4]-1)*self._seq[i].radius/1000.0,1))+" km altitude)")
        print(
            "Beta plane angle:        " + str(round(x[6 + (i - 1) * 4]*180/pi,0))+"° ")
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
        

    print("\nArrival at " + self._seq[-1].name)
    print("Arrival epoch:           " + str(round(t_P[-1].mjd2000*Edy2Kdy,1)) + " KUT ("+Kdate(t_P[-1].mjd2000*Edy2Kdy,planet_pack)+")")
    print("Total mission time:      " + str(round(sum(T)*Edy2Kdy,1)) + " days ("+Ktime(sum(T)*Edy2Kdy,planet_pack)+")")   
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
    if self._orbit_insertion:
        # In this case we compute the insertion DV as a single pericenter
        # burn. We assume that the final orbit is coplanar with the incoming hyperbola plane.
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
        
    return {
        "objective_dv": DV_Sum,
        "ejection_dv": ejection_dv,
        "ejection_strategy": ejection["strategy"],
        "ejection_direct_dv": ejection["direct_dv"],
        "ejection_split_dv": ejection["split_dv"],
        "ejection_soi_correction_dv": ejection["soi_correction_dv"],
        "ejection_dv_full": DV_Ej[0],
        "ejection_dv_xy": DV_Ej[1],
        "ejection_dv_z": DV_Ej[2],
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


