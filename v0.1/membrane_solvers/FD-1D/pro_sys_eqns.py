# -*- coding: utf-8 -*-
"""
Created on Tue Jul 24 21:33:32 2018

Water FLux Equation

@author: Drizdar
"""
import math
import operator
from scipy.optimize import brentq
import pandas as pd
import numpy as np
import os


#%%##Water Flux Equation#######################################################
'''
Jw = flux across membrane for PRO power generation
iterative equation, solve using brentq function (fzero function in MATLAB)
F should be equal to 0
'''
def wfe(J_w,k,D,C_D,C_F,S,B,A,dP,C_Vant,T,p_D,p_F):
    T = T + 273.15 #Converts T to Kelvin
    F = lambda J_w: J_w -(A*(C_Vant*(T)*((p_D*C_D*math.exp(-(J_w/k))- \
    p_F*C_F*math.exp((J_w*S*0.001)/D))/(1+(B/J_w)*(math.exp((J_w*S*0.001)/D) \
    -math.exp(-(J_w/k)))))-dP))
    J_w = brentq(F,0.0001,300)
    return J_w

#%%#Saltwater Density Equation#################################################
'''
From "Millero, F., Poisson, A., 1981. International one-atmosphere equation of state of seawater. Deep Sea Res Part Oceanogr Res Pap 28, 625–629."
Equation vaid from 0 - 40C, and from 0.5 to 43 ppt, standard error is 3.6e-3 kg/m^3
note, will be using this equation for brine until high salinity density function can be found
S should be in g/kg or ppt
T should be in C
'''
def swden(S,T):
    #S = float(S)
    A = 0.824493 - 0.0040899*T + 0.000076438*T**2 - 0.00000082467*T**3 + 0.0000000053875*T**4
    B = -0.00572466 + 0.00010227*T - 0.0000016546*T**2
    C = 0.00048314
    #p_o is the density of water, from "Bigg, 1967, British Journal of Applied Physics, 8, 521-537"
    p_o = 999.842594 + 0.06793952*T - 0.00909529*T**2 + 0.0001001685*T**3 - 0.000001120083*T**4 + 0.000000006536332*T**5 #kg/m^3
    p_SW = p_o + A*S + B*S**1.5 + C*S**2 #kg/m^3
    return p_SW

#%%##Membrane Foul Equation####################################################
def mfe(mu, A_max):
    #R_m = membrane specific resistance
    R_m = (1/(A_max*mu))
    return R_m

#%%##Water Perm Equation#######################################################
def wpe(mu, R_m, R_f):
    #A = water permeability across the membrane
    AF = (1/mu)*(1/(R_m + R_f))
    return AF

#%%#Salt Flux Equation#########################################################
def sfe(J_w,k,D,C_D,C_F,S,B,A):
    #J_s = salt flux across the membrane
    #K = S/D
    #J_s = ((B*0.001)/(K*(B*0.001)+1)*(C_D-C_F-K*J_w*0.001*C_F))*1000 #g L m-2 h-1 kg-1 - alternate eqn form
    J_s = B*((C_D*math.exp(-(J_w/k)) - C_F*math.exp((J_w*S*0.001)/D))\
             /(1+(B/J_w)*(math.exp((J_w*S*0.001)/D) - math.exp(-(J_w/k)))))
    return J_s

#%%#Salinity Converter#########################################################
'''
Converts salinity from mol/L to g/kg (ppt), T in celsius
not really used, but good to have
'''
def salcon(M,T,MW_NaCl):
    F = lambda s: (s - ((M * MW_NaCl * 1000)/(swden(s,T))))
    sal = brentq(F,0.001,1000) 
    return sal

#%%#Colebrook Equation#########################################################
def Colebrook(ks,D,rho,V,mu):
    f = 0
    '''
    ks = td.ks_vals['Drawn tubing']
    D = 0.004 #m
    rho = 1.23 #kg/m^3
    V = 50 #m/s
    mu = 1.79e-5 #N-s/m^2
    '''
    Re = (rho * V * D)/mu
    if Re < 2000:
        f = 64/Re
    else:
        F = lambda f: (1/math.sqrt(f)) + 2 * math.log10((ks/(3.7*D*1000)) + (2.51/(Re*math.sqrt(f))))
        f = brentq(F,0.00001,100)
    return f

#%%#Energy Equation############################################################
def engeq(ks,D,rho,mu,z1,z2,L,g,Q,Pump,inf,kLtot):
    'used to find amount of power consumed by and cost of pumping for a given segment'
    'Assumes that velocity is equal at beginning and end'
    'Assume 0 pressure at beginning, system needs to be pressurized to 20 psi to maintain system pressure,'
    'as stated on pg 47 in "EPRI, 2013. Electricity Use and Management in the Municipal Water Supply and Wastewater Industries ( No. 3002001433). Electric Power Research Institute."'
    pr2 = 1.38e5 #Pa ~ 20 psi
    A = (D**2)*(math.pi/4) #m^2
    #Q = Q * FD_r #Assumes transmitted water is used as feed solution
    V = Q/A #m/s
    f = Colebrook(ks,D,rho,V,mu)
    kLtot = sum(kLtot)
    phd = (z2 - z1) + (f*(L/D)+kLtot)*((V**2)/(2*g)) + (pr2/(rho*g)) #m
    P = ((rho*g*Q*phd)/Pump.W_n)/1000 #kW
    Q = Q * 3600 #m^3/hr
    Tr_use = P/Q #kWh/m^3
    Cpl = (362.34 * (D**2) + 78.5*D + 15.478)*1.2273 #$/m, from Texas Cost Doc, adjusted from Mar 2007 - Nov 2018
    Pipe_Cost = L * Cpl #$
    if P > 0:
        Pump_Cost = (Pump.a * (P*1.341) ** Pump.n) * inf #1.341 converts kW to horsepower
    else:
        Pump_Cost = 0
    Cost = Pipe_Cost + Pump_Cost #$
    return Tr_use, phd, P, D, V, Pipe_Cost, Pump_Cost, Cost, Pump.a, Pump.n
    #for labor, assume a pay rate
    
    #note: in the EPRI reference, this eqn was used for pumping energy cost:
   # Electricity (kWh/day) = ((Flow(gpm) * pumping head (in feet))/(3960 * pumping efficiency)) * 0.746 * 24
   #average pumping efficiency of 65% was assumed

#%%#Hazen-Williams Equation####################################################
#Not currently used, but here in case needed
def HazenWilliams(L,Q,D,C):
    #From "Hammer, M.J., Hammer, Jr., M.J., 2012. Water and Wastewater Technology, 7th ed. Pearson Education Inc., Upper Saddle River, NJ."
    hL = 0.002131 * L*((100/C)**1.85) * ((Q**1.85)/(D**4.8655)) #m
    return hL

#%%#Optimal Water Flux and Energy Density######################################   
def owfed(k,D,C_D,C_F,S,B,A,C_Vant,T,p_D,p_F):
    dP_vec  = [] #bar
    W_vec   = [] #Watt/m^2
    J_w_vec = [] #L m-2 h-1
    n   = 10000
    
    for i in range(0,n):
        dP_vec.append((i+1)/10)
        dPi = dP_vec[i] #bar
        J_w_vec.append(wfe(J_w_vec,k,D,C_D,C_F,S,B,A,dPi,C_Vant,T,p_D,p_F)) #L m-2 h-1
        w = J_w_vec[i] * dP_vec[i] * (1/36) #Watt/m^2
        W_vec.append(w) #Watt/m^2
        if J_w_vec[i] <= 0.3:
            break
   
    index, W_max = max(enumerate(W_vec), key=operator.itemgetter(1)) #finds maximum value and index in a list
    maxW_dP = dP_vec[index]
    maxW_J_w = J_w_vec[index]
    J_w = maxW_J_w #L m-2 h-1 bar-1
    W = W_max #Watt m-2
    dP = maxW_dP #bar
    return J_w, W, dP

#%%#Financial Projection#######################################################
def finproj(N,EP,dEP,r,MW_net,OpTime,Cost,Cost_OM):
    'calculates overall net present value and ROI for a PRO facility'
    #N is # of years, EP is energy price ($/kw), dEP is Energy price increase,
    #MW_net is net Power (MW), Optime is operation time (%)
    #Vrp is present value of revenue
    #need to figure out realistic operation time
    #for now, just do 75% of the year
    thryr = 24*365*OpTime #hr/yr - operation time for plant
    E_net = thryr * MW_net #MWh/yr
    Vrp = 0 #Present value of revenue
    for n in range(1,int(N+1)):
        Vrpn = (E_net * 1000 * EP * (1+dEP)**n) * (1/((1+r)**n)) - Cost_OM
        Vrp = Vrp + Vrpn
        if Vrp >= Cost:
            PB_Pd = n
            break
    Vrp = 0
    for n in range(1,int(N+1)):
        Vrpn = (E_net * 1000 * EP * (1+dEP)**n) * (1/((1+r)**n)) - Cost_OM
        Vrp = Vrp + Vrpn
    if Cost >= Vrp:
        PB_Pd = -404 #Means that payback period is either greater than plant lifetime, or that Vrp is negative
    return Vrp,PB_Pd,E_net
'''
    from desalting handbook:
    interest rate of 5%, service life of 25 years
    indirect costs as follows (% of construction costs):
        freight and insurance: 5%
        construction overhead and profit: 15%
        owners cost: 10%
        contingency cost: 10%
        membrane replacement: 14% - assumes average membrane life of >5 years
'''
#%%#Specific Energy Calculator(Ideal)##########################################
#Based on Equations from "Lin, S., Straub, A., Elimelech, M., 2014. Thermodynamic limits of extractable energy by pressure retarded osmosis. Energy Environ Sci 7, 2706–2714."
def SEcalc(CD,CF,v,T,R,kind):
    T = T + 273.15 #convert temp to Kelvin
    if kind == "rev": #Thermodynamically reversible process
        SE = math.exp(((CD*math.log(CD)-CF*math.log(CF))/(CD-CF)) - 1) - \
        ((CD*CF)/(CD - CF)) * (math.log(CD) - math.log(CF))
    elif kind == "CoCur": #Co-current circulatory process
        SE = ((math.sqrt(CD) - math.sqrt(CF))**2)/4
    elif kind == "CntCur": #Counter-current circulatory process
        SE = 0.25 *((CD-CF)**2/(CD+CF))
    SE = SE * v * T * R
    return SE

#%%#PT Energy Calculator#######################################################
def pte(PT_d,PT_f,PTopp,Q_f,Q_d,Q_m):
    #PTopp: 0 = no pretreatment, 1 = Draw Pretreatment, 2 = Feed pretreatment, 3 = feed and draw pretreatment
    if PTopp == 0:
        SE_PT = 0
        P_PT = 0
        #Q_PT = 0
        ptcd = 0
        ptcf = 0
    elif PTopp == 1:
        SE_PT = PT_d.PT_E * (Q_d/Q_m) 
        P_PT = (PT_d.PT_E * Q_d)/1e6 #MW
        #multiplies by ratio so that energy consumption from pt is only for portion of water that needs pt
        #Q_PT = Q_d
        ptcd = 1
        ptcf = 0
    elif PTopp == 2:
        SE_PT = PT_f.PT_E *  (Q_f/Q_m)
        P_PT = (PT_f.PT_E * Q_f)/1e6 #MW
        #Q_PT = Q_f
        ptcd = 0
        ptcf = 1
    elif PTopp == 3:
        SE_PT = (PT_d.PT_E * (Q_d/Q_m)) + (PT_f.PT_E *  (Q_f/Q_m))
        P_PT = ((PT_d.PT_E * Q_d)/1e6) + ((PT_f.PT_E * Q_f)/1e6) #MW
        #Q_PT = Q_m
        ptcd = 1
        ptcf = 1
    return SE_PT, ptcd, ptcf, P_PT

#%%#Membrane System Performance Evaluator######################################
def comboin(MemTyp,TurbTyp,Tr_usei,Q,C_D,C_F,proj,OpTime,Cms,T,mu,Tropp,PT_d,PT_f,C_Vant,PTopp,inf,v,R,M_geometry,MW_NaCl):
    import Mem_Flux
    'note: memTyp is from the props.MembraneType class, TurbType is from the props.TurbineType class'
    #sys_in = owfed(MemTyp.k,MemTyp.D,C_D,C_F,MemTyp.S,MemTyp.B,MemTyp.A,C_Vant)
    sys_in = Mem_Flux.CntCurMod(C_Vant,MemTyp.k,MemTyp.D,MemTyp.S,MemTyp.B,MemTyp.A,M_geometry,C_D,C_F,T,mu)
    #Tuple contains  Q_f[end],Q_d[end],Q_f[start],Q_d[start],W_avg,Se_avg,kw_gross,C_f[end],C_d[end],J_w_avg,dP[start],J_s_avg,SysCon,dP[end]
    #gross > net
    A_m = M_geometry[6]
    Mpv = M_geometry[8]
    Memb_LT = 7 #years, average membrane lifetime - assume complete system replacement every X years
    #V_mod = 0.03142 #m^3 Volume of an individual module, from "Baker, R.W., 2004. Membrane Technology and Applications, 2nd ed. John Wiley and Sons, Ltd."
    #Pd_mod = 775 #m^2/m^3 Packing Density, from "Sundaramoorthy, S., Srinivasan, G., Murthy, D.V.R., 2011. An analytical model for spiral wound reverse osmosis membrane modules: Part I — Model development and parameter estimation. Desalination 280, 403–411."
    #MWt_NaCl = 58.442769 #g/mol
    Q = Q*3600*1000 #L/hr - Q draw influent total
    #A_t = (Q / sys_in[3])*A_m*Mpv #m^2
    #V_pv = V_mod * 8 #Volume of a pressure vessel
    n_protr = (Q / sys_in[3])/14 #number of process trains (skids)
    n_protr = math.ceil(n_protr) #round up to find total number of process trains (skids)
    A_t = n_protr*14*Mpv*A_m #calculate total membrane area
    n_mods = A_t/(A_m) #number of membrane modules/elements
    n_pv = n_mods/Mpv #number of pressure vessels, logic from applied membranes.com (for 8in dia membrane, pressure vessels hold max of 7 elements)
    Q_di = sys_in[3] * n_pv #L/hr - total draw influent
    Q_fi = sys_in[2] * n_pv #L/hr
    Q_mi = Q_di + Q_fi #L/hr total inflow rate
    Q_de = sys_in[1] * n_pv #L/hr - total effluent rate from draw 
    dQ_d = Q_de - Q_di #effluent sent to turbine - remainder is sent to PX to maintain system pressure
    Q_fe = sys_in[0] * n_pv #L/hr total unused feed solution (try to minimize this)
    Tr_usei.Q = Q_fi/(3600*1000) # converts transmission flow into feed influent in m3/s
    Tr_use = Tr_usei.calc()
    aa = Tr_use[8] #extract pump cost coefficients
    nn = Tr_use[9]
    Tr_use = Tr_use * Tropp #calculates energy used in transmission and multiplies values by Tropp 
    #Tr_use output is Tr_use, phd, P, D, V, Pipe_Cost, Pump_Cost, Cost, Pump.a, Pump.n
    pt = pte(PT_d,PT_f,PTopp,Q_fi,Q_di,Q_mi)
    #pte output is SE_PT, ptcd, ptcf, P_PT
    SE_PT = pt[0] #Energy Used by Pretreatment 
    C_d = (C_D * swden(C_D,T)) / MW_NaCl #mol/m^3
    C_f = (C_F * swden(C_F,T)) / MW_NaCl #mol/m^3
    C_de = sys_in[8] #ppt (mg/kg)
    C_fe = sys_in[7] #ppt (mg/kg)
    C_me = (C_fe*Q_fe + C_de*Q_de)/(Q_fe+Q_de) #concentration of mixed effluent
    SE_max_rev = SEcalc(C_d,C_f,v,T,R,"rev") #kWh/m^3, Maximum Specific Energy for a Thermodynamically reversible process
    SE_max_CoCur = SEcalc(C_d,C_f,v,T,R,"CoCur") #kWh/m^3, Maximum Specific Energy for a Co-current circulatory process
    SE_max_CntCur = SEcalc(C_d,C_f,v,T,R,"CntCur") #kWh/m^3, Maximum Specific Energy for a Counter-current circulatory process
    #SE_gross2 = sys_in[1]/sys_in[0] #kWh/m^3, W/J_w - Q_f based instead of Q_m based
    SE_gross = sys_in[5] #kWh/m^3
    #SE_gross3 = (sys_in[2]*Q_f)/(Q_m*36) #kWh/m^3
    Circ = (((sys_in[10]-sys_in[13]*0.97)*Q_di)*(1/36)*(1e-6))/0.89 #MW - energy used to repressurize draw solution 
    SE_net = SE_gross*TurbTyp.W_n - Tr_use[0]*(Q_fi/Q_mi) - SE_PT - Circ/Q_mi #kWh/m^3
    #MW_gross2 = (A_t*sys_in[1])/1e6 #MW #only works of Q_d = Q_f
    MW_gross = (dQ_d * sys_in[13]) *(1/36) *(1e-6) #dQ_d * dP[end]
    #MW_gross2 = (SE_gross*Q_de)/1e6 - old way
    #MW_net2 = (SE_net*Q_de)/1e6 #MW - old way
    #Circ2 = (1.4 * Q_mi * (1/36)*(1e-6))/0.9 #energy used for circulating the two streams at system pressure of 20 psi - note: assume not needed since water already pressurized by external systems coming in
    MW_net = MW_gross*TurbTyp.W_n - (Tr_use[2]/1000) - pt[3] - Circ
    if MW_net < 0:
        MW_net = 0 #to make the math work, since MW_net < 0 gives an imaginary number for Turbcost
    Cost_Memb = A_t *Cms # cost in dollars
    #Yearly operation/maintainence cost - based on membrane replacement schedule - assume half cost since Pressure vessels do not need to be replaced (only membrane elements)
    Cost_OM = (Cost_Memb/Memb_LT)/2 #$/Yr
    Cost_Tr = Tr_use[7] #$
    #turb cost MW multiplied by 2 to account for cost of both turbine and PX
    #turb cost multiplied by 1.1865 to account for inflation from November 2008 to November 2018. Inflation calc from https://www.bls.gov/data/inflation_calculator.htm
    #Turbine cost multiplied by 1.26 to account for exchange rate between Euros and $ in November 2008. Sourced from: https://freecurrencyrates.com/en/exchange-rate-history/EUR-USD/2008/cbr
    Cost_Turb = (2600 * (( MW_net * 1000)**0.54) * 1.26) *1.1865 #$ eqn 19b from "Aggidis, G.A., Luchinskaya, E., Rothschild, R., Howard, D.C., 2010. The costs of small-scale hydro power production: Impact on the development of existing potential. Renew Energ 35, 2632–2638."
    Cost_Turb = Cost_Turb + (2600 * (( Circ * 1000)**0.54) * 1.26) *1.1865 #adds PX
    Cost_Turb = Cost_Turb + (aa * (Circ*1.341*1000) ** nn) * inf #1.341 * 1000 converts MW to horsepower - adds booster pump
    #for pt: PT_E, ptcd, ptc
    Cost_PTd = pt[1]*PT_d.costQ(Q_di,inf) #$
    Cost_PTf = pt[2]*PT_f.costQ(Q_fi,inf) #$
    Cost_PT = Cost_PTd + Cost_PTf #$
    Cost_Cons = (Cost_Memb + Cost_Tr + Cost_Turb + Cost_PT)*1.35 #multiply by 1.35 to add cost of pipes (0.1), controls (0.2), and site prep (0.05), as suggested by McGiveney
    Cost = Cost_Cons * 1.35 #multiply by 1.35 to add Engineering, Legal, and Administrative Costs, as suggested by McGiveney
    if SE_net <=0 or MW_net <= 0:
        #placeholder values since trial would have a negative MW_net, which is needed for PV calculations
        Cost_Turb = -404
        Cost_Unit = -404 
        PV_rev = -404
        PB_pd = -404
        PV_net = -404
        E_net = -404
        Cost_OM_Unit = -404 
        Cost_OM_Unit2 = -404 
        SysCon2 = 0 #Indicates that trial would have a negative SE_net
    else:
        Cost_Unit = Cost/(MW_net * 1000)  #$/kW, compare to EIA Cost Estimates
        PV_revv = finproj(proj[0],proj[1],proj[2],proj[3],MW_net,OpTime,Cost,Cost_OM) #$, gross savings of implementation
        PB_pd = PV_revv[1]
        PV_rev = PV_revv[0]
        PV_net = PV_rev - Cost #$, net energy cost savings, aka, the ROI
        E_net = float(PV_revv[2]) #MWh/Yr
        Cost_OM_Unit = Cost_OM/(MW_net*1000) #$/kw-yr, Compare to EIA Cost Estimates
        Cost_OM_Unit2 = Cost_OM/(E_net*1000) #$/kWh, compare to Megaton Project
        SysCon2 = 1 #Indicates trial would have a positive SE_net
    SE_Ineff  = SE_gross * (1-TurbTyp.W_n) + Circ/Q_mi
    SE_Tr = Tr_use[0]*(Q_fi/Q_mi)
    MW_Ineff = MW_gross*(1 - TurbTyp.W_n) + Circ
    MW_PT = (SE_PT*Q_de)/1e6
    MW_Tr = (Tr_use[0]*Q_de)/1e6
    SysCon = sys_in[12]
    #A_t = (Q_feed/C_n) / J_w
    #n_m = A_t/A_m
    #SE_gross = W/J_w - units = kWh/m^3 (triple checked, math checks out)
    #SE_net = SE_gross*W_n-Tr_use-PT_use
    #MW_gross = (A_t*W)/1e6
    #MW_net = (SE_net*Q_feed)/1e6
    #Ineff_E = SE_gross * (1-W_n)
    #Ineff_W = (Ineff_E*Q_feed)/1e6
    #PT_W = (PT_use*Q_feed)/1e6
    return sys_in[9],sys_in[11],sys_in[4],sys_in[10],sys_in[13],A_t,n_pv,n_protr,C_de,C_fe,C_me,Q_di,Q_fi,Q_de,Q_fe,SE_max_rev,SE_max_CoCur,SE_max_CntCur,SE_gross,SE_net,MW_gross,Circ,MW_net,E_net,\
            Cost_Memb,Cost_Tr,Cost_Turb,Cost_PT,Cost_Cons,Cost,Cost_Unit,Cost_OM,Cost_OM_Unit,Cost_OM_Unit2,PV_rev,PV_net,PB_pd,SE_Ineff,SE_PT,SE_Tr,MW_Ineff,MW_PT,MW_Tr,SysCon, SysCon2

def combo(membs,TurbTyp,Tr_usei,Q,C_D,C_F,proj,OpTime,Cms,T,mu,Tropp,PT_d,PT_f,C_Vant,PTopp,inf,v,R,M_geometry,MW_NaCl):
    p = pd.Categorical(['J_w(L m-2 h-1)','J_s(g L m-2 h-1 kg-1)','W(W m-2)','dPi(bar)','dPe(bar)','A_t(m^2)','n_pv','n_protr','C_de(ppt)','C_fe(ppt)','C_me(ppt)','Q_di(L/hr)','Q_fi(L/hr)','Q_de(L/hr)','Q_fe(L/hr)',\
            'SE_max_rev(kWh/m^3)','SE_max_CoCur(kWh/m^3)','SE_max_CntCur(kWh/m^3)','SE_gross(kWh/m^3)','SE_net(kWh/m^3)',\
            'MW_gross(MW)','Circ','MW_net(MW)','E_net(MWh/yr)','Cost_Memb($)','Cost_Tr($)','Cost_Turb($)','Cost_PT($)','Cost_Cons($)','Cost($)',\
            'Cost_Unit($/kW)','Cost_OM($/Yr)','Cost_OM_Unit($/kw-yr)','Cost_OM_Unit2($/kwh)','PV_rev($)','PV_net($)','PB_pd(Yrs)','SE_Ineff(kWh/m^3)','SE_PT(kWh/m^3)','SE_Tr(kWh/m^3)','MW_Ineff(MW)','MW_PT(MW)','MW_Tr(MW)','SysCon','SysCon2'])
    sys_out = pd.DataFrame()
    it = 1
    for i in range(0,len(membs)):
        print('it = %s'%it)
        j = comboin(membs[i],TurbTyp,Tr_usei,Q,C_D,C_F,proj,OpTime,Cms,T,mu,Tropp,PT_d,PT_f,C_Vant,PTopp,inf,v,R,M_geometry,MW_NaCl)
        it = it + 1
        f = membs[i].name
        sys_out['%s'%f] = pd.Series(j)
    #sys_out['Units'] = pd.Categorical(['L m-2 h-1','W m-2','bar','m^2','#','kWh/m^3','kWh/m^3', 'MW','MW','$'])
    sys_out.index = p
    sys_out = sys_out.transpose()
    return sys_out   

#%%#Comparison Case Generator##################################################
def Compcg(sw_ww,roc_ww,roc_sw,varcomp,mmx):
    #can find data pertaining to min and max values for a given parameter
    #Filters out trials with negative net power
    sw_ww = sw_ww[sw_ww.SysCon2 == 1]
    roc_ww = roc_ww[roc_ww.SysCon2 == 1]
    roc_sw = roc_sw[roc_sw.SysCon2 == 1]
    
    if mmx == 'max':
        roc_sw_var = roc_sw.loc[roc_sw[varcomp].idxmax()]
        roc_ww_var = roc_ww.loc[roc_ww[varcomp].idxmax()]
        sw_ww_var = sw_ww.loc[sw_ww[varcomp].idxmax()]
        membid = pd.DataFrame({"Membrane":[roc_sw[varcomp].idxmax(),roc_ww[varcomp].idxmax(),sw_ww[varcomp].idxmax()]},
                                        index = ["sw_ww","roc_ww","roc_sw"])
    elif mmx == 'min':
        roc_sw_var = roc_sw.loc[roc_sw[varcomp].idxmin()]
        roc_ww_var = roc_ww.loc[roc_ww[varcomp].idxmin()]
        sw_ww_var = sw_ww.loc[sw_ww[varcomp].idxmin()]
        membid = pd.DataFrame({"Membrane":[roc_sw[varcomp].idxmin(),roc_ww[varcomp].idxmin(),sw_ww[varcomp].idxmin()]},
                                        index = ["sw_ww","roc_ww","roc_sw"])

    comp_var1 = pd.DataFrame(data=[sw_ww_var,
                              roc_ww_var,
                              roc_sw_var]
                        ,index=["sw_ww","roc_ww","roc_sw"])
    comp_var = membid.join(comp_var1)
    return comp_var

#%%#sensitivity and uncertainty analysis#######################################
'''Returns all of the info as both a .csv file and a Pandas Dataframe'''
def comboUSA(prob,N_Yrs,Tropp,PT_d,PT_f,C_Vant,PTopp,inf,v,R,M_geometry,g,MW_NaCl,csv_name):
    import TabularData as td
    import props
    #membs = prob[:,0:5]
    k = prob[:,0]
    D = prob[:,1]
    S = prob[:,2]
    B = prob[:,3]
    A = prob[:,4]
    Turbeff = prob[:,5]
    z1 = prob[:,6]
    z2 = prob[:,7]
    L = prob[:,8]
    d = prob[:,9]
    Q = prob[:,10]
    C_D = prob[:,11]
    C_F = prob[:,12]
    proj = np.copy(prob[:,12:16]) #do np.copy to prevent overwrite in next line
    proj[:,0] = N_Yrs #Sets the first value in Proj to N, since N is Fixed
    OpTime = prob[:,16]
    Cms = prob[:,17]
    T = prob[:,18]
    
    Pmp = td.Pump_dat
    
    #P is same as P for normal analysis --> whenever edit is made to other P, edit this one too
    p = pd.Categorical(['J_w(L m-2 h-1)','J_s(g L m-2 h-1 kg-1)','W(W m-2)','dPi(bar)','dPe(bar)','A_t(m^2)','n_pv','n_protr','C_de(ppt)','C_fe(ppt)','C_me(ppt)','Q_di(L/hr)','Q_fi(L/hr)','Q_de(L/hr)','Q_fe(L/hr)',\
            'SE_max_rev(kWh/m^3)','SE_max_CoCur(kWh/m^3)','SE_max_CntCur(kWh/m^3)','SE_gross(kWh/m^3)','SE_net(kWh/m^3)',\
            'MW_gross(MW)','Circ','MW_net(MW)','E_net(MWh/yr)','Cost_Memb($)','Cost_Tr($)','Cost_Turb($)','Cost_PT($)','Cost_Cons($)','Cost($)',\
            'Cost_Unit($/kW)','Cost_OM($/Yr)','Cost_OM_Unit($/kw-yr)','Cost_OM_Unit2($/kwh)','PV_rev($)','PV_net($)','PB_pd(Yrs)','SE_Ineff(kWh/m^3)','SE_PT(kWh/m^3)','SE_Tr(kWh/m^3)','MW_Ineff(MW)','MW_PT(MW)','MW_Tr(MW)','SysCon','SysCon2'])
    sys_out = pd.DataFrame()
    it = 1
    l = len(prob)
    for i in range(0,len(Turbeff)): #len Turbeff is arbitraty, any length will work
        print('it = {0}/{1}'.format(it,l))
        #Create a synthetic membrane
        memb = props.MembraneType('%s'%it,'%s'%it,k[i],D[i],S[i],B[i],A[i],"material","ref")
        #Create a synthetic turbine for analysis
        TurbTyp = props.TurbineType("Turb", Turbeff[i])
        #Create a synthetic pipe network - # bends, material, and pipe size kept same as Tampa case for consistency
        # eqn below (rho) from -> from "Jones, F.E., Harris, G.L., 1992. ITS-90 density of water formulation for volumetric standards calibration. Journal of Research of the National Institute of Standards and Technology 97."
        rho = 999.84847 + 6.337563e-2 * T[i] - 8.523829e-3 * (T[i]**2) + 6.943248e-5 * (T[i]**3) - 3.821216e-7 * (T[i]**4) #kg/m^3, valid from 5 to 40 C, for air-saturated fresh water
        mu  = 0.001*(1 + 0.636*(T[i]-20)/41)**(-1/0.636)   #Kg m-1 s-1 -> from "Pawlowski, J., 1991. Veränderliche Stoffgrössen in der Ähnlichkeitstheorie. Frankfurt am Main: Salle."
        Tr_usei = props.Tr_use("Desal Plant","Wastewater Plant",
                      td.ks_vals['PVC and Plastic Pipes'],rho,d[i],mu,z1[i],z2[i],L[i],g,Q[i],Pmp["Sub"],inf,
                      [td.kL_vals['Entrance']['Slightly rounded'],td.kL_vals['Exit']['Slightly rounded'],
                       td.kL_vals['Elbows']['Long radius 90°, flanged']*4,td.kL_vals['Elbows']['Long radius 45°, flanged']*2])
        j = comboin(memb,TurbTyp,Tr_usei,Q[i],C_D[i],C_F[i],proj[i],OpTime[i],Cms[i],T[i],mu,\
                    Tropp,PT_d,PT_f,C_Vant,PTopp,inf,v,R,M_geometry,MW_NaCl)
        it = it + 1
        f = memb.name
        sys_out['%s'%f] = pd.Series(j)
    sys_out.index = p
    sys_out = sys_out.transpose()
    path = 'SA tests' #folder that test files will be saved to
    sys_out.to_csv(os.path.join(path,r'%s.csv'%csv_name))
    return sys_out 

#%%#Sensitivity Analysis Input Saver###########################################
'''Saves inputs from Sensitivity Analysis as dataframe and CSV File'''
def SAin(problem,param_values,csv_name):
    
    p = pd.Categorical(['k(L m-2 h-1)', 'D(m^2/h)', 'S(m)', 'B(L m-2 h-1)', 'A(L m-2 h-1 bar-1)', 'Turbeff',\
    'z1(m)', 'z2(m)', 'L(m)', 'd(m)', 'Q(m^3/s)', 'C_D(g/kg)', 'C_F(g/kg)', 'EP($/kWh)', 'dEP($)', 'r', 'OpTime', 'Cms($/m^2)','T'])
    
    sys_in = pd.DataFrame(data=param_values,
              index=np.arange(1, (len(param_values)+1)),
              columns=p)
    inbound = pd.DataFrame.from_dict(problem['bounds'])
    if len(p) != len(problem['names']) == True:
        inbound.index = pd.DataFrame.from_dict(problem['names'])
    else:
        inbound.index = p
    inbound.columns = ['min','max']
    inbound = inbound.transpose()

    path = 'SA tests' #folder that test files will be saved to
    inbound.to_csv(os.path.join(path,r'%s.csv'%csv_name))
    
    sys_in.to_csv(os.path.join(path,r'%s.csv'%csv_name),mode ='a')
    return sys_in

#%%#Sensitivity Analysis Output Saver##########################################
'''Saves inputs from Sensitivity Analysis as CSV File'''
#Type: input the analysis type, e.g. "sobol"

def SA_indices(problem,Si,Y,method,csv_title_SI,Notes,N_Samples):
    path = 'SA tests' #folder that test files will be saved to
    M = pd.DataFrame({"Method used: %s"%method,"Sensitivity Analysis on %s"%Y,"Notes: %s"%Notes,
                      "Number of samples: %s"%N_Samples})
    M.to_csv(os.path.join(path,r'%s.csv'%csv_title_SI),mode ='a',header=False,index=False)
    
    if method.lower() in ["rbd"]:
        DS1 = pd.DataFrame({"","S1"})
        DS1.to_csv(os.path.join(path,r'%s.csv'%csv_title_SI),mode ='a',header=False,index=False)
        S1 = pd.DataFrame.from_dict(Si['S1'])
        S1.index = problem['names']
        S1.to_csv(os.path.join(path,r'%s.csv'%csv_title_SI),mode ='a',header=False)
    
    if method.lower() in ["fast"]:
        DS1 = pd.DataFrame({"","S1"})
        DS1.to_csv(os.path.join(path,r'%s.csv'%csv_title_SI),mode ='a',header=False,index=False)
        S1 = pd.DataFrame.from_dict(Si['S1'])
        S1.index = problem['names']
        S1.to_csv(os.path.join(path,r'%s.csv'%csv_title_SI),mode ='a',header=False)
        
        DST = pd.DataFrame({"","ST"})
        DST.to_csv(os.path.join(path,r'%s.csv'%csv_title_SI),mode ='a',header=False,index=False)
        ST = pd.DataFrame.from_dict(Si['ST'])
        ST.index = problem['names']
        ST.to_csv(os.path.join(path,r'%s.csv'%csv_title_SI),mode ='a',header=False)
        
    elif method.lower() in ["sobol"]:
        DS1 = pd.DataFrame({"","S1"})
        DS1.to_csv(os.path.join(path,r'%s.csv'%csv_title_SI),mode ='a',header=False,index=False)
        S1 = pd.DataFrame.from_dict(Si['S1'])
        S1.index = problem['names']
        S1.to_csv(os.path.join(path,r'%s.csv'%csv_title_SI),mode ='a',header=False)
        
        DS1_conf = pd.DataFrame({"","S1_conf"})
        DS1_conf.to_csv(os.path.join(path,r'%s.csv'%csv_title_SI),mode ='a',header=False,index=False)
        S1_conf = pd.DataFrame.from_dict(Si['S1_conf'])
        S1_conf.index = problem['names']
        S1_conf.to_csv(os.path.join(path,r'%s.csv'%csv_title_SI),mode ='a',header=False)
        
        DS2 = pd.DataFrame({"","S2"})
        DS2.to_csv(os.path.join(path,r'%s.csv'%csv_title_SI),mode ='a',header=False,index=False)
        S2 = pd.DataFrame.from_dict(Si['S2'])
        S2.index = problem['names']
        S2.columns = problem['names']
        S2.to_csv(os.path.join(path,r'%s.csv'%csv_title_SI),mode ='a')
        
        DS2_conf = pd.DataFrame({"","S2_conf"})
        DS2_conf.to_csv(os.path.join(path,r'%s.csv'%csv_title_SI),mode ='a',header=False,index=False)
        S2_conf = pd.DataFrame.from_dict(Si['S2_conf'])
        S2_conf.index = problem['names']
        S2_conf.columns = problem['names']
        S2_conf.to_csv(os.path.join(path,r'%s.csv'%csv_title_SI),mode ='a')
        
        DST = pd.DataFrame({"","ST"})
        DST.to_csv(os.path.join(path,r'%s.csv'%csv_title_SI),mode ='a',header=False,index=False)
        ST = pd.DataFrame.from_dict(Si['ST'])
        ST.index = problem['names']
        ST.to_csv(os.path.join(path,r'%s.csv'%csv_title_SI),mode ='a',header=False)
        
        DST_conf = pd.DataFrame({"","ST_conf"})
        DST_conf.to_csv(os.path.join(path,r'%s.csv'%csv_title_SI),mode ='a',header=False,index=False)
        ST_conf = pd.DataFrame.from_dict(Si['ST_conf'])
        ST_conf.index = problem['names']
        ST_conf.to_csv(os.path.join(path,r'%s.csv'%csv_title_SI),mode ='a',header=False)
        
#%%#Param Values Generator#####################################################
'''Used for individual SA'''
def PVG(param_value,B_Vals,switch):
    param_values = np.zeros(shape = (len(param_value),19))
    param_values[:,0] = B_Vals['k'][switch]
    param_values[:,1] = B_Vals['D'][switch]
    param_values[:,2] = B_Vals['S'][switch]
    param_values[:,3] = B_Vals['B'][switch]
    param_values[:,4] = B_Vals['A'][switch]
    param_values[:,5] = B_Vals['Turbeff'][switch]
    param_values[:,6] = B_Vals['z1'][switch]
    param_values[:,7] = B_Vals['z2'][switch]
    param_values[:,8] = B_Vals['L'][switch]
    param_values[:,9] = B_Vals['d'][switch]
    param_values[:,10] = B_Vals['Q'][switch]
    param_values[:,11] = B_Vals['C_D'][switch]
    param_values[:,12] = B_Vals['C_F'][switch]
    param_values[:,13] = B_Vals['EP'][switch]
    param_values[:,14] = B_Vals['dEP'][switch]
    param_values[:,15] = B_Vals['r'][switch]
    param_values[:,16] = B_Vals['OpTime'][switch]
    param_values[:,17] = B_Vals['Cms'][switch]
    param_values[:,18] = B_Vals['T'][switch]
    return param_values

