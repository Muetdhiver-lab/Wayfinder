# -*- coding: utf-8 -*-
"""
Created on Wed May 29 08:02:50 2019

@author: v.fave
"""

"""
This is a monkey patch of pykep's mga_1dsm 
The functions are dedicated to decoding and displaying
the mission plan (in terminal)
"""

from numpy import dot,cross,linalg,around,rad2deg, array
from pykep import epoch
from scipy.linalg import norm
from pykep.core import DAY2SEC, lambert_problem, propagate_lagrangian, fb_prop, ic2par
from math import sqrt, pi, cos, sin, acos, log, atan, asin, copysign
from _Vanilla_System import body as Vanilla_body,row as Vanilla_row,col as Vanilla_col
from _JNSQ_System import body as JNSQ_body,row as JNSQ_row,col as JNSQ_col
import numpy as np

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
        print("VSOI[1] is too small, set to 1e-12")
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
    theta = acos((a * (1 - e * e) - Rsoi) / (e * Rsoi))
    return theta

def EJ_angle_to_Prograde(self,Pe,Pro):

    pro_zeroZ = np.array([Pro[0],Pro[1],0])
    pro_zeroZ = pro_zeroZ/linalg.norm(pro_zeroZ)
    if cross(Pe,pro_zeroZ)[2] < 0 :
        return 2*pi - acos(dot(Pe,pro_zeroZ))
    else :
        return acos(dot(Pe,pro_zeroZ))

def EJ_details(self,v_EJDV_at_SOI,Rsoi,prograde_uv,Alt):
    EJDV_at_SOI = linalg.norm(v_EJDV_at_SOI)
    theta = self.EJ_angle_from_Pe(EJDV_at_SOI,Rsoi,Alt)
    #print(str(theta*180/pi)+" EJ angle from Pe")
    ejectionDirection = v_EJDV_at_SOI/ EJDV_at_SOI

    if (abs(sin(theta)) < abs(ejectionDirection[2])) :
        print("Warning, theta : "+str(round(theta,1))+" and ejdr[2] : "+str(round(ejectionDirection[2],1)))
        ejectionDeltaV = 99999,99999,99999
        return [ejectionDeltaV, float(90), float(90)]
    else :
        periapsisDirection = self.EJ_Pe_Direction(ejectionDirection, theta)
        ejectionAngle = self.EJ_angle_to_Prograde(periapsisDirection, prograde_uv)
        Pe_X_Ej = cross(periapsisDirection, ejectionDirection)
        ejectionInclination = acos((Pe_X_Ej/linalg.norm(Pe_X_Ej))[2])
        ejectionInclination = copysign(ejectionInclination,pi - theta)
        ejectionInclination = copysign(ejectionInclination,ejectionDirection[2])
        ejectionDeltaV = self.EJ_BurnDV(v_EJDV_at_SOI,Alt,Rsoi,ejectionInclination,Alt)
        return ejectionDeltaV, ejectionInclination, ejectionAngle;
'''
def deco_fitness_function(self, x):
        
    theta = 2 * pi * x[1]
    phi = acos(2 * x[2] - 1) - pi / 2
    Vinfx = x[3] * cos(phi) * cos(theta)
    Vinfy = x[3] * cos(phi) * sin(theta)
    Vinfz = x[3] * sin(phi)
    Vinfxy = sqrt(Vinfy**2+Vinfx**2)
        
    FirstBody = self._seq[0]
    Alt = 100000        
    Rsoi = FirstBody.orbital_elements[0] * (FirstBody.mu_self/FirstBody.mu_central_body)**(0.4)
    v0 = sqrt(FirstBody.mu_self / (FirstBody.radius+Alt) )
    print(v0)
    vxy_ob = sqrt(Vinfxy**2 + 2 * FirstBody.mu_self / (FirstBody.radius+Alt) -2 * FirstBody.mu_self / Rsoi) - v0 # oberth only for xy
    print("vxy_ob = "+str(vxy_ob))
    v_ob = sqrt(vxy_ob**2+Vinfz**2)
    print("v_ob = "+str(v_ob))
    print("v_esc_100 = "+str(sqrt(FirstBody.mu_self / FirstBody.radius+Alt)))
    print("v_esc_RSoi = "+str(sqrt(FirstBody.mu_self / Rsoi)))

    if v_ob < (vxy_ob+abs(Vinfz)):
        deco_fitness = self.fitness(x)[0]+v_ob-x[3]
    else :
        deco_fitness = self.fitness(x)[0]+vxy_ob+Vinfz-x[3]        
    print(self.fitness(x)[0]+vxy_ob+Vinfz-x[3] )        
    #print("orig fitness : "+str(orig_fitness_function(self,dv)))
    #print("new  fitness : "+str(fitness))        

    return deco_fitness
'''


def transx(self,x,alt=100000,planet_pack="Vanilla"):
    if planet_pack == "Vanilla":
        body = Vanilla_body
        row  = Vanilla_row
        col  = Vanilla_col
        Edy2Kdy = 4
    elif planet_pack == "JNSQ":
        body = JNSQ_body
        row  = JNSQ_row
        col  = JNSQ_col
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
    Rsoi = body[row[self._seq[0].name],col['R_soi (km)']]*1000

    DV_Ej,DV_i = self.EJ_details(Vinf,Rsoi,fward_P[0],alt)[0:2]
    DV_i = DV_i * 180/pi
    #4 - depending on case, split the ej burn or not.
    if DV_Ej[0] < (DV_Ej[1]+Vinfz):
        DV_Sum = DV_Ej[0]
    else :
        DV_Sum = DV_Ej[1]+abs(Vinfz)
        
    pe = alt+self._seq[0].radius
    a = -1*self._seq[0].mu_self/(x[3]*x[3])
    e = (a-pe)/a

    Eta = acos(-1/e)
    Eta_Rsoi = acos((a * (1 - e * e) - Rsoi) / (e * Rsoi));
    Delta = asin(1/e)*2
    Gamma = atan(dot(r_P[0],Vinf)/(norm(Vinf)*norm(r_P[0])))
    
    v0 = [a + b for a, b in zip(v_P[0], [Vinfx, Vinfy, Vinfz])]
    
    k_el = ic2par(r_P[0],v0,self.common_mu) 
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

    r, v = propagate_lagrangian(
        r_P[0], v0, x[4] * T[0] * DAY2SEC, self.common_mu)

    print("DSM after :              " + Ktime(x[4] * T[0]*Edy2Kdy,planet_pack))

    # Lambert arc to reach seq[1]
    dt = (1 - x[4]) * T[0] * DAY2SEC
    l = lambert_problem(
        r, r_P[1], dt, self.common_mu, cw = False, max_revs=self.max_revs)
    v_end_l = l.get_v2()[0]
    v_beg_l = l.get_v1()[0]

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
        v_out = fb_prop(v_end_l, v_P[i], x[
                7 + (i - 1) * 4] * self._seq[i].radius, x[6 + (i - 1) * 4], self._seq[i].mu_self)
        print(
            "Fly-by epoch:            " + Kdate(t_P[i].mjd2000*Edy2Kdy,planet_pack)+" ("+Ktime(t_P[i].mjd2000*Edy2Kdy-t_P[0].mjd2000*Edy2Kdy,planet_pack)+" after launch)")
        print(
            "Fly-by radius:           " + str(round(x[7 + (i - 1) * 4],2)) + " planetary radii ("+str(round((x[7 + (i - 1) * 4]-1)*self._seq[i].radius/1000.0,1))+" km altitude)")
        print(
            "Beta plane angle:        " + str(round(x[6 + (i - 1) * 4]*180/pi,0))+"° ")
        # s/c propagation before the DSM
        r, v = propagate_lagrangian(
            r_P[i], v_out, x[8 + (i - 1) * 4] * T[i] * DAY2SEC, self.common_mu)
        print("DSM after :              " + Ktime(x[8 + (i - 1) * 4] * T[i]*Edy2Kdy,planet_pack) + " after flyby ("+Ktime((x[8 + (i - 1) * 4] * T[i]+t_P[i].mjd2000-t_P[0].mjd2000)*Edy2Kdy,planet_pack)+" days after launch)")
        # Lambert arc to reach Earth during (1-nu2)*T2 (second segment)
        dt = (1 - x[8 + (i - 1) * 4]) * T[i] * DAY2SEC
        l = lambert_problem(r, r_P[i + 1], dt, self.common_mu, cw=False, max_revs=self.max_revs)
        v_end_l = l.get_v2()[0]
        v_beg_l = l.get_v1()[0]
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
        
        
def decode_dV_tof(self,x,alt=100000,planet_pack="Vanilla"):
    if planet_pack == "Vanilla":
        body = Vanilla_body
        row  = Vanilla_row
        col  = Vanilla_col
        Edy2Kdy = 4
    elif planet_pack == "JNSQ":
        body = JNSQ_body
        row  = JNSQ_row
        col  = JNSQ_col
        Edy2Kdy = 2
        
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
    DV_Ej = self.EJ_details(Vinf,Rsoi,fward_P[0],alt)[0]

    if DV_Ej[0] < (DV_Ej[1]+Vinfz):
        DV_Sum = DV_Ej[0]
    else :
        DV_Sum = DV_Ej[1]+abs(Vinfz)

    v0 = [a + b for a, b in zip(v_P[0], [Vinfx, Vinfy, Vinfz])]
    r, v = propagate_lagrangian(
        r_P[0], v0, x[4] * T[0] * DAY2SEC, self.common_mu)


    '''Lambert arc to reach seq[1]'''
    dt = (1 - x[4]) * T[0] * DAY2SEC
    l = lambert_problem(r, r_P[1], dt, self.common_mu, cw = False, max_revs=self.max_revs)
    v_end_l = l.get_v2()[0]
    v_beg_l = l.get_v1()[0]

    '''First DSM occuring at time nu1*T1'''
    DV[0] = norm([a - b for a, b in zip(v_beg_l, v)])

    DV_Sum += DV[0]
    
    '''And we proceed with each successive leg'''
    for i in range(1, self.n_legs):
        # Fly-by
        v_out = fb_prop(v_end_l, v_P[i], x[7 + (i-1) * 4] * self._seq[i].radius, x[6 + (i-1) * 4], self._seq[i].mu_self)
        # s/c propagation before the DSM
        r, v = propagate_lagrangian(
            r_P[i], v_out, x[8 + (i-1) * 4] * T[i] * DAY2SEC, self.common_mu)
        # Lambert arc to reach Earth during (1-nu2)*T2 (second segment)
        dt = (1 - x[8 + (i-1) * 4]) * T[i] * DAY2SEC
        l = lambert_problem(r, r_P[i + 1], dt, self.common_mu, cw=False, max_revs=self.max_revs)
        v_end_l = l.get_v2()[0]
        v_beg_l = l.get_v1()[0]
        # DSM occuring at time nu2*T2
        DV[i] = norm([a - b for a, b in zip(v_beg_l, v)])
        DV_Sum += DV[i]   

    DV[-1] = norm([a - b for a, b in zip(v_end_l, v_P[-1])])
    if self._orbit_insertion:
        # In this case we compute the insertion DV as a single pericenter
        # burn. We assume that the final orbit is coplanar with the incoming hyperbola plane.
        DVper = np.sqrt(DV[-1] * DV[-1] + 2 *
                                 self._seq[-1].mu_self / self._rp_target)
        DVper2 = np.sqrt(2 * self._seq[-1].mu_self / self._rp_target -
                                 self._seq[-1].mu_self / self._rp_target * (1. - self._e_target))
        DV_Inj = np.abs(DVper - DVper2)    
            
    if self._orbit_insertion:
        DV_Sum += DV_Inj
    elif self._add_vinf_arr  :
        DV_Sum += DV[-1]
    return DV_Sum, round(t_P[0].mjd2000*Edy2Kdy,1) ,round(sum(T)*Edy2Kdy,1), norm(Vinf)


