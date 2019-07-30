# -*- coding: utf-8 -*-
"""
Created on Mon Oct 15 12:03:22 2018

Membrane Flux Model

Basis: For a given membrane size, calculate ideal Q for system using co-current 
approximation, and then use given system size for counter-current approximation. 

@author: Drizdar
"""

import numpy as np
import math
import pro_sys_eqns as pse

#def Mem_Flux(MemTyp,FD_r,C_D,C_F,C_Vant,A_m):
#Calculates dP, dQF, dQD Jw,Js,CD,CF
#This is for a co-current flow situation
#Based on equations from "Straub, A., Lin, S., Elimelech, M., 2014. Module-scale analysis of pressure retarded osmosis: performance limitations and implications for full-scale operation. Environ Sci Technology 48, 12435–44."

#Parameters
"""
v = 2 #Van't Hoff factor for NaCl, from "Lin, S., Straub, A., Elimelech, M., 2014. Thermodynamic limits of extractable energy by pressure retarded osmosis. Energy Environ Sci 7, 2706–2714."
R = 0.083144598 #L-bar/mol-K
MW_NaCl = 58.442769 #g/mol
C_Vant = (v*R/MW_NaCl) #L-bar/g-K - Van't Hoff Coefficient (Modified)

k   = 99 #L m-2 h-1, mass transfer coefficient
D   = 5.328e-6 #m^2/h, diffusion constant
S   = 0.00031 #m, structural parameter
B   = 0.09 #L m-2 h-1, salt permeability of the active layer
A   = 5.11 #L m-2 h-1 bar-1, water permeability of the active layer

A_m = 4.18 #m^2
L_m = 1.016 #m - membrane element length (A)
Mpv = 1 #number of membrane elements per pressure vessel
C_D = 28 #g/kg, solute concentration on draw side
C_F = 0.01 #g/kg, solute concentration on feed side
T = 22 #deg C
"""

def CoCurMod(C_Vant,k,D,S,B,A,A_m,L_m,Mpv,C_D,C_F,T):
    p_D = pse.swden(C_D,T)/1000 #kg/L
    p_F = pse.swden(C_F,T)/1000 #kg/L
    Ld = A_m / L_m #m^2/m, represents linear density of membrane surface -> multiply Q values by this amount to get total flow capacity
    #Number of elements
    N = 301
    A_mt = A_m * Mpv #multiply area by number of elements per pressure vessel
    L_mt = L_m * Mpv #multiply length by number of elements per pressure vessel
    s_vec = np.linspace(0,L_mt,N) #start at 0, increment to A_m and the length of the vector is N
    ds = s_vec[2]-s_vec[1] #m^2
    
    #Discretize
    C_d = np.zeros((N,1)) #a 1 by N array
    C_f = np.zeros((N,1))
    Q_d = np.zeros((N,1))
    Q_f = np.zeros((N,1))
    J_w = np.zeros((N,1))
    J_s = np.zeros((N,1))
    #W = np.zeros((N,1))
    #Se = np.zeros((N,1))
    #Mw = np.zeros((N,1))
    
    #Initialize
    C_d[0] = C_D #g/kg Starting condition
    C_f[0] = C_F #value at end - use to check - for now, do [0,0], next, move on to [0,N-1]
    Q_d[0] = 12000/Ld #L/hr
    Q_f[0] = 12000/Ld #L/hr
    
    sys_in = pse.owfed(k,D,C_d[0],C_f[0],S,B,A,C_Vant,T,p_D,p_F)
    dP = sys_in[2] #bar
    
    def CoCur(Q_d,Q_f,C_d,C_f,J_w,J_s,N,dP,p_D,p_F):
        for n in range(0,N-1):
            J_w[n] = pse.wfe(J_w[n],k,D,C_d[n],C_f[n],S,B,A,dP,C_Vant,T,p_D,p_F)
            J_s[n] = pse.sfe(J_w[n],k,D,C_d[n],C_f[n],S,B,A)
            
            Q_d[n+1] = J_w[n]*ds + Q_d[n]
            Q_f[n+1] = -J_w[n]*ds + Q_f[n]
            
            C_d[n+1] = (-J_s[n]*ds + Q_d[n]*C_d[n])/Q_d[n+1]
            C_f[n+1] = (J_s[n]*ds + Q_f[n]*C_f[n])/Q_f[n+1]
            
            p_D = pse.swden(C_d[n+1],T)/1000 #kg/L
            p_F = pse.swden(C_f[n+1],T)/1000 #kg/L
            
            
            #W[n+1] = (dP*(Q_d[n]-Q_d[0])/s_vec[n])/36
            #Se[n+1] = ((dP*(Q_d[n]-Q_d[0]))/(Q_d[0]+Q_f[0]))/36 #kWh/m^3
            #Mw[n+1] = (Q_d[n] * Se[n])/1e3 #in kW
            #if J_s[0,n] >= J_w[0,n]:
                #break
        J_w[N-1] = pse.wfe(J_w[N-1],k,D,C_d[N-1],C_f[N-1],S,B,A,dP,C_Vant,T,p_D,p_F)
        J_s[N-1] = pse.sfe(J_w[N-1],k,D,C_d[N-1],C_f[N-1],S,B,A)
    
    CoCur(Q_d,Q_f,C_d,C_f,J_w,J_s,N,dP,p_D,p_F)
    
    #Adjust starting flows so that none of the Q_f is wasted
    it = 1
    #print('it = %s'%it)
    
    while Q_f[N-1] >= 0:
        Q_f[0] = Q_f[0] - (Q_f[N-1]/10)
        Q_d[0] = Q_f[0]
        CoCur(Q_d,Q_f,C_d,C_f,J_w,J_s,N,dP,p_D,p_F)
        it = it + 1
        #print('it = %s'%it)
        if J_s[0] >= J_w[0]:
            break
        if math.isclose(J_w[N-1],(0.05*J_w[0]),rel_tol = 0.25) == True or J_s[N-1] > 5*J_w[N-1]:
                break
            
    Q_f = np.multiply(Q_f,Ld) #scale Q by linear density
    Q_d = np.multiply(Q_d,Ld)
    
    W_avg = float((dP*(Q_d[N-1]-Q_d[0])/A_mt)/36) #W/m^2 - for system
    J_w_avg = np.mean(J_w)
    J_s_avg = np.mean(J_s)
    
    #W3 = sys_in[1] - same as W2, just different form
    
    Se_avg = float(((dP*(Q_d[N-1]-Q_d[0]))/(Q_d[0]+Q_f[0]))/36) #kWh/m^3
    kw_gross = float((Q_d[N-1] * Se_avg)/1e3)
    """
    W2 = float(J_w[0] * dP * (1/36)) #W/m^2 - max based on coupon scale system
    print("CoCur system actual values: SE_avg = %f kWh/m^3,"%Se_avg,"kW_gross = %f kW,"%kw_gross,"W_avg = %f W/m^2"%W_avg )
    #SE_gross = (W2*A_mt)/(Q_d[0]+Q_f[0]) #Q_m based - not really valid since A_mt won't have same W
    Se_gross_max = float(((dP*Q_f[0])/(Q_d[0]+Q_f[0]))/36) #max for 100% transfer
    #SE_gross3 = (W_avg*A_mt)/(Q_d[0]+Q_f[0]) #same as Se_avg
    kw_gross_max = float((Q_d[N-1] * Se_gross_max)/1e3)
    print("Cocur system max theoretical values: SE_max %f kWh/m^3,"%Se_gross_max,"kW_gross_max = %f kW,"%kw_gross_max,"W_coupon = %f W/m^2 (for coupon scale membrane)"%W2 )
    print("Cocur % of max: SE =",((Se_avg/Se_gross_max)*100),"%, kW_gross = ",((kw_gross/kw_gross_max)*100),"%, W_avg = ",((W_avg/W2)*100),"%")
    """
    return int(np.round(Q_f[0])),int(np.round(Q_d[0])),int(np.round(Q_f[N-1])),int(np.round(Q_d[N-1])),float(J_s[0]),float(J_w[0]),W_avg,Se_avg,kw_gross,float(C_f[N-1]),float(C_d[N-1]),float(J_w_avg),dP,float(J_s_avg)


def CntCurMod(C_Vant,k,D,S,B,A,A_m,L_m,Mpv,C_D,C_F,T):
    CoCur = CoCurMod(C_Vant,k,D,S,B,A,A_m,L_m,Mpv,C_D,C_F,T)
    Ld = A_m / L_m #m^2/m, represents linear density of membrane surface -> multiply Q values by this amount to get total flow capacity
    if CoCur[4] > CoCur[5]:
        SysCon = 0
        #Tuple contains  Q_f[end],Q_d[end],Q_f[start],Q_d[start],W_avg,Se_avg,kw_gross,C_f[end],C_d[end],J_w_avg,dP,SysCon
        return CoCur[2],CoCur[3],CoCur[0],CoCur[1],CoCur[6],CoCur[7],CoCur[8],CoCur[9],CoCur[10],CoCur[11],CoCur[12],CoCur[13],SysCon
    
    N = 301
    A_mt = A_m * Mpv #multiply area by number of elements per pressure vessel
    L_mt = L_m * Mpv #multiply length by number of elements per pressure vessel
    s_vec = np.linspace(0,L_mt,N) #start at 0, increment to A_m and the length of the vector is N
    ds = s_vec[2]-s_vec[1] #m^2
    
    Ccc = 0 #CoCurrent Check  - if 1, then return Cocurrent numbers
    #Discretize
    p_D = pse.swden(C_D,T)/1000 #kg/L
    p_F = pse.swden(C_F,T)/1000 #kg/L
    
    C_d = np.zeros((N,1)) #a 1 by N array
    C_f = np.zeros((N,1))
    Q_d = np.zeros((N,1))
    Q_f = np.zeros((N,1))
    J_w = np.zeros((N,1))
    J_s = np.zeros((N,1))
    #W = np.zeros((N,1))
    #Se = np.zeros((N,1))
    #Mw = np.zeros((N,1))
    
    #Initialize
    C_fchk = C_F #final C_f value should be this
    Q_fchk = (2*CoCur[0])/Ld #final Q_f value should be this - set to beginning Qf value from cocurrent
    """
    C_fchk = C_F #final C_f value should be this
    Q_fchk = 12000 #final Q_f value should be this
    C_d[0] = C_D #g/kg Starting condition
    C_f[0] = 2.8
    Q_d[0] = 12000 #L/hr
    Q_f[0] = 294 #L/hr
    """
    C_d[0] = C_D #g/kg Starting condition
    C_f[0] = CoCur[9] #set to final CF value from cocurrent
    Q_d[0] = 2*CoCur[1]/Ld #L/hr set to initial QD value from cocurrent
    Q_f[0] = CoCur[2]/Ld #L/hr set to final QF value from cocurrent
    
    #sys_in = pse.owfed(k,D,C_d[0],C_fchk,S,B,A,C_Vant,T,p_D,p_F)
    #dP = sys_in[2] #bar
    dP = CoCur[12] #bar
    p_F = pse.swden(C_f[0],T)/1000 #kg/L #for starting feed density
    
    def CntCur(Q_d,Q_f,C_d,C_f,J_w,J_s,N,dP,Ccc,p_D,p_F):
        for n in range(0,N-1):
            J_w[n] = pse.wfe(J_w[n],k,D,C_d[n],C_f[n],S,B,A,dP,C_Vant,T,p_D,p_F)
            J_s[n] = pse.sfe(J_w[n],k,D,C_d[n],C_f[n],S,B,A)
            
            Q_d[n+1] = J_w[n]*ds + Q_d[n]
            Q_f[n+1] = J_w[n]*ds + Q_f[n]
            
            C_d[n+1] = (-J_s[n]*ds + Q_d[n]*C_d[n])/Q_d[n+1]
            C_f[n+1] = (-J_s[n]*ds + Q_f[n]*C_f[n])/Q_f[n+1]
            
            p_D = pse.swden(C_d[n+1],T)/1000 #kg/L
            p_F = pse.swden(C_f[n+1],T)/1000 #kg/L

            if C_f[n+1] < -0.05:
                Ccc = 1
                return Q_d,Q_f,C_d,C_f,J_w,J_s,Ccc #means that countercurrent won't work (given current numerics)
            
            #W[n+1] = (dP*(Q_d[n]-Q_d[0])/s_vec[n])/36
            #Se[n+1] = ((dP*(Q_d[n]-Q_d[0]))/(Q_d[0]+Q_f[0]))/36 #kWh/m^3
            #Mw[n+1] = (Q_d[n] * Se[n])/1e3 #in kW
            
        J_w[N-1] = pse.wfe(J_w[N-1],k,D,C_d[N-1],C_f[N-1],S,B,A,dP,C_Vant,T,p_D,p_F)
        J_s[N-1] = pse.sfe(J_w[N-1],k,D,C_d[N-1],C_f[N-1],S,B,A)
        
        return Q_d,Q_f,C_d,C_f,J_w,J_s,Ccc
        
    CntCur(Q_d,Q_f,C_d,C_f,J_w,J_s,N,dP,Ccc,p_D,p_F)
    
    if Ccc == 1:
        SysCon = 0
        #Tuple contains  Q_f[end],Q_d[end],Q_f[start],Q_d[start],W_avg,Se_avg,kw_gross,C_f[end],C_d[end],J_w_avg,dP,SysCon,J_s_avg
        return CoCur[2],CoCur[3],CoCur[0],CoCur[1],CoCur[6],CoCur[7],CoCur[8],CoCur[9],CoCur[10],CoCur[11],CoCur[12],CoCur[13],SysCon
            
    it = 1
    #it2 = 1
    #print('it2 = %s'%it2)
      
    #print("almost done")
    while math.isclose(C_f[N-1],C_fchk,rel_tol = 0.1) == False\
            and math.isclose(Q_f[N-1],Q_fchk,rel_tol = 0.5) == False:
        C_f[0] = C_f[0] + ((C_fchk - C_f[N-1])/2)
        p_F = pse.swden(C_f[0],T)/1000 #kg/L #for starting feed density
        Q_f[0] = Q_f[0] + ((Q_fchk - Q_f[N-1])/10)
        CntCur(Q_d,Q_f,C_d,C_f,J_w,J_s,N,dP,Ccc,p_D,p_F)
        it = it + 1
        #print('it = %s'%it)
        if math.isclose(C_f[N-1],C_fchk,rel_tol = 0.1) == True\
            and math.isclose(Q_f[N-1],Q_fchk,rel_tol = 0.5) == True:
            it = 0
            break
    
    """
    while Q_f[0] >= 0:
        Q_fchk = Q_f[N-1] - Q_f[0]
        Q_d[0] = Q_fchk
        print('it2 = %s'%it2)
        it2 = it2 + 1
        while C_f[N-1] != C_fchk and Q_f[N-1] != Q_fchk:
            C_f[0] = C_f[0] + ((C_fchk - C_f[N-1])/2)
            Q_f[0] = Q_f[0] + ((Q_fchk - Q_f[N-1])/10)
            CntCur(Q_d,Q_f,C_d,C_f,J_w,J_s,N,dP,Ccc)
            it = it + 1
            print('it = %s'%it)
            if (math.isclose(C_f[N-1],C_fchk,rel_tol = 0.1) == True\
                and math.isclose(Q_f[N-1],Q_fchk,rel_tol = 1) == True) or\
                math.isclose(J_w[N-1],0.5,rel_tol = 0.20) == True or J_s[N-1] > 5*J_w[N-1]\
                or math.isclose((Q_f[0]/Q_f[N-1]),0.01,rel_tol = 1) == True: 
                it = 0
                break
        if math.isclose(J_w[N-1],0.5,rel_tol = 0.20) == True or J_s[N-1] > 5*J_w[N-1]\
            or math.isclose((Q_f[0]/Q_f[N-1]),0.01,rel_tol = 1) == True:
                break
    """
        
    #print('done')
    
    Q_f = np.multiply(Q_f,Ld) #scale Q by linear density
    Q_d = np.multiply(Q_d,Ld)
    
    
    W_avg = float((dP*(Q_d[N-1]-Q_d[0])/A_mt)/36) #W/m^2 - for system
    #W3 = sys_in[1] - same as W2, just different form
    J_w_avg = np.mean(J_w)
    J_s_avg = np.mean(J_s)
    
    Se_avg = float(((dP*(Q_d[N-1]-Q_d[0]))/(Q_d[0]+Q_f[N-1]))/36) #kWh/m^3
    kw_gross = float((Q_d[N-1] * Se_avg)/1e3)
    """
    W2 = float(J_w[0] * dP * (1/36)) #W/m^2 - max based on coupon scale system
    #W3 = sys_in[1] - same as W2, just different form
    kw_gross2 = (W_avg * A_mt)/1000 #should be in Kw - makes more sense - use this one - more representative of gross
    print("system actual values: SE_avg = %f kWh/m^3,"%Se_avg,"kW_gross = %f kW,"%kw_gross2,"W_avg = %f W/m^2"%W_avg )
    #SE_gross = (W2*A_mt)/(Q_d[0]+Q_f[0]) #Q_m based - not really valid since A_mt won't have same W
    Se_gross_max = float(((dP*Q_f[N-1])/(Q_d[0]+Q_f[N-1]))/36) #max for 100% transfer
    Se_gross_max2 = float(((dP)/36)) #max for 100% transfer
    #SE_gross3 = (W_avg*A_mt)/(Q_d[0]+Q_f[N-1]) #same as Se_avg
    kw_gross_max = float((Q_d[N-1] * Se_gross_max)/1e3)
    #kw_gross_max2 = (W2*A_mt)/1000
    print("system max theoretical values: SE_max %f kWh/m^3,"%Se_gross_max,"kW_gross_max = %f kW,"%kw_gross_max,"W_coupon = %f W/m^2 (for coupon scale membrane)"%W2 )
    print("% of max: SE =",(Se_avg/Se_gross_max),"%, kW_gross = ",(kw_gross/kw_gross_max),"%, W_avg = ",(W_avg/W2),"%")
    """
    SysCon = 1
    return int(np.round(Q_f[0])),int(np.round(Q_d[N-1])),int(np.round(Q_f[N-1])),int(np.round(Q_d[0])),W_avg,Se_avg,kw_gross,float(C_f[0]),float(C_d[N-1]),J_w_avg,dP,J_s_avg,SysCon

#TestThingn = CntCurMod(C_Vant,k,D,S,B,A,A_m,L_m,Mpv,C_D,C_F,T) #Copy and paste this in if you want it to run - for some reason, doesn't work otherwise