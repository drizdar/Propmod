# -*- coding: utf-8 -*-
"""
Created on Tue Jul 24 21:23:17 2018

Basic Plant Designer / Process Model

@author: Drizdar

Desalination energy recovery system designer
"""

import pro_sys_eqns as pse
import props
import TabularData as td
import pandas as pd
import plotter
from timeit import default_timer as timer
     

#%%#Universal Input Parameters#################################################
#It would be neat to make each of these a dictionary with the different recorded
#Values, as well as where they came from. - These are constants - dictionary not needed
start = timer() #used for timing the function

#C_Vant = 0.7345 #bar kg g-1 -> Van't Hoff Coefficient, from He, W., Wang, Y., Shaheed, M., 2015. Stand-alone seawater RO (reverse osmosis) desalination powered by PV (photovoltaic) and PRO (pressure retarded osmosis). Energy 86, 423–435.

v = 2 #Van't Hoff factor for NaCl, from "Lin, S., Straub, A., Elimelech, M., 2014. Thermodynamic limits of extractable energy by pressure retarded osmosis. Energy Environ Sci 7, 2706–2714."
R = 0.083144598 #L-bar/mol-K
MW_NaCl = 58.442769 #g/mol
C_Vant = (v*R/MW_NaCl) #L-bar/g-K - Van't Hoff Coefficient (Modified)
R = 2.30957E-06 #kWh/Mol-K #Change Units of R
#Note: Inflation value for Turbine (Found in pse) may also need to be adjusted depending on when model is ran
A_m = 40.88 #m^2 ~440 sq. ft, from VI email chain
L_m = 1.016 #m - membrane element length (A)
Mpv = 7 #number of membrane elements per pressure vessel
h_c = 0.0007 #m - channel height
l_f = 0.0045 #m - length between spacers

T = 31.6 #deg C
mu  = 0.001*(1 + 0.636*(T-20)/41)**(-1/0.636)   #Kg m-1 s-1 -> from "Pawlowski, J., 1991. Veränderliche Stoffgrössen in der Ähnlichkeitstheorie. Frankfurt am Main: Salle."
# eqn below (rho) from -> from "Jones, F.E., Harris, G.L., 1992. ITS-90 density of water formulation for volumetric standards calibration. Journal of Research of the National Institute of Standards and Technology 97."
rho = 999.84847 + 6.337563e-2 * T - 8.523829e-3 * (T**2) + 6.943248e-5 * (T**3) - 3.821216e-7 * (T**4) #kg/m^3, valid from 5 to 40 C, for air-saturated fresh water
#Cms = 14.18 #Converted from Au $ to USD under exchange rate of 0.709 USD/Au $ From "Helfer, F., Lemckert, C., 2015. The power of salinity gradients: An Australian example. Renewable and Sustainable Energy Reviews 50, 1–16."
Cms = 5.80 #$/m^2, from "Naghiloo, A., Abbaspour, M., Mohammadi-Ivatloo, B., Bakhtari, K., 2015. Modeling and design of a 25 MW osmotic power plant (PRO) on Bahmanshir River of Iran. Renew Energ 78, 51–59."
#Cms = 12.23 #$/m2 from VI email chain (440 ft^2 per element, around $500 each - email says 400 ft2, but Dow website says 440, so I'll go with Dow)
g = 9.81 #m/s^2
inf = 1.2089 #inflation from September 2007 to November 2018 for cost estimates, from https://www.bls.gov/data/inflation_calculator.htm


#%%#Membrane and spacer Types#############################################################
membs = td.membs

memdat = pd.DataFrame.from_records([s.to_dict() for s in membs]) #DataFrame with all membrane data - not in order
memdat.set_index("0_Name", inplace=True)

membs[0].displayStats() 

sd = td.Spacer_dat
Spacers = pd.DataFrame.from_records([sd[s].to_dict(h_c,l_f) for s, value in sd.items()]) #values in to_dict are h_c, l_f, in m
Spacers.set_index("Shape", inplace=True)

shape = "1:3 ellipse"
M_geometry = [Spacers.loc[shape,"alpha"],Spacers.loc[shape,"beta"],1,Spacers.loc[shape,"d_h"],
              h_c, l_f, A_m, L_m, Mpv] #values are alpha, beta, gamma, d_h, h_c, l_f, A_m, L_m, Mpv
          
#%%#Plants and Objects#########################################################
#Create Plants
#Desal plants contain C_ROC (high), C_SW (low), Q_ROC(m^3/s),Optime,Elev
#WW plants contain C_Eff (low), C_SW (high), Q_Eff (m^3/s),Optime,Elev
#note: Tampa Bay Salinity is assumed to be 27 ppt
TBdp = props.DesalPlant("Tampa Bay Desalination Plant",66,27,0.83,0.80,3.1) #66,27,0.83,0.6,3.1
HCww = props.WWPlant("Howard F. Curren Wastewater Plant",0.431,27,2.29,1,2.7)
SCRA = props.WWPlant("South County Regional AWWTP",0.256,27,0.2,1,13.4)

#for seven seas water plants, missing elevation data
SSWdp = props.DesalPlant("Seven Seas Water Desalination Plant",66,27,0.83,0.75,3.1)
SSWww = props.WWPlant("Seven Seas Water Wastewater Plant",0.01,27,2.29,1,2.7)

#Create Turbines
Pelton = props.TurbineType("Pelton Turbine", 0.85)
Hydro = props.TurbineType("Hydro-Turbine", 0.90)
RPE = props.TurbineType("Rotary Pressure Exchanger", 0.97)
#Pump = props.PumpType("Pump", 0.75) #Value from "Hammer, M.J., Hammer, Jr., M.J., 2012. Water and Wastewater Technology, 7th ed. Pearson Education Inc., Upper Saddle River, NJ."

Pmp = td.Pump_dat
Pmp_dat = pd.DataFrame.from_records([Pmp[s].to_dict() for s, value in Pmp.items()]) #DataFrame with all Pump data - not in order
Pmp_dat.set_index("ID", inplace=True)

#Calculate Pretreatment Energy Consumption in kWh/m^3 
Pt = td.Pt_dat

Pt_dat = pd.DataFrame.from_records([Pt[s].to_dict() for s, value in Pt.items()]) #DataFrame with all Pretreatment data - not in order
Pt_dat.set_index("ID", inplace=True)

#Calculate Transmission Energy Consumption and Cost
#note: props.Tr_use contains (DesalPlant,WWPlant,ks,rho(kg/m^3),D(m),mu(kg/m-s),z1(m),z2(m),L(m),g(m/s^2),Q(l/hr),Pump,inf,kLtot)
#Note: establish base conditions for input - Q will change depending on membrane configuration
HCww_TBdpi = props.Tr_use("Tampa Bay Desalination Plant","Howard F. Curren Wastewater Plant",
                      td.ks_vals['PVC and Plastic Pipes'],rho,1,mu,HCww.Elev,TBdp.Elev,23500,g,TBdp.Q_ROC,Pmp["Sub"],inf,
                      [td.kL_vals['Entrance']['Slightly rounded'],td.kL_vals['Exit']['Slightly rounded'],
                       td.kL_vals['Elbows']['Long radius 90°, flanged']*5,td.kL_vals['Elbows']['Long radius 45°, flanged']*6])
SCRA_TBdpi = props.Tr_use("Tampa Bay Desalination Plant","South County Regional AWWTP",
                      td.ks_vals['PVC and Plastic Pipes'],rho,0.46,mu,SCRA.Elev,TBdp.Elev,9382,g,SCRA.Q_Eff,Pmp["Sub"],inf,
                      [td.kL_vals['Entrance']['Slightly rounded'],td.kL_vals['Exit']['Slightly rounded'],
                       td.kL_vals['Elbows']['Long radius 90°, flanged']*4,td.kL_vals['Elbows']['Long radius 45°, flanged']*2])
#HCww_TBdp = HCww_TBdpi.calc()
#print()
#HCww_TBdpi.displayStats()
#print(HCww_TBdp)

#proj sets financial calc parameters - N (Max # years), EP (Energy price, $/kwh),
#dEP (change in Energy Price per year, $), r (interest rate)
#see the payback period - will be acceptable if less than 7 years
#Set maximum years for payback period to 30
proj = [50,0.09,0.003,0.03]

#%%#Plant Comparer#############################################################
'''
note: format is pse.combo(membs,TurbTyp,Tr_usei,Q,C_D,C_F,proj,OpTime,Cms,T,mu\
Tropp,PT_d,PT_f,C_Vant,PTopp,inf,v,R,M_geometry,MW_NaCl
#PTopp: 0 - no pretreatment, 1 - Draw Pretreatment, 2 = Feed pretreatment, 3 = feed and draw pretreatment
#Tropp = 1 = transmission, 0 = no transmission
#Assume ROC and WW do not need pretreatment - place Pt['MF'] in as a default though to make code work
#note: multiply Tr_use by 0 if no transmission necessary
#PV_rev and PV_net represent present value of revenue and net present value (revenue - cost) over given time-frame (first # in proj())
'''
#ROC-WW using the HCww as the WW Plant
#roc_ww = pse.combo(membs,RPE,HCww_TBdpi,TBdp.Q_ROC,TBdp.C_ROC,HCww.C_Eff,proj,TBdp.OpTime,Cms,T,mu,\
#                   1,Pt['MF'],Pt['MF'],C_Vant,0,inf,v,R,M_geometry,MW_NaCl)

#ROC-WW trial using SCRA as the WW Plant
roc_ww2 = pse.combo(membs,RPE,SCRA_TBdpi,SCRA.Q_Eff*1.4,TBdp.C_ROC,SCRA.C_Eff,proj,TBdp.OpTime,Cms,T,mu,\
                   1,Pt['MF'],Pt['MF'],C_Vant,0,inf,v,R,M_geometry,MW_NaCl)

sw_ww = pse.combo(membs,Hydro,HCww_TBdpi,HCww.Q_Eff,HCww.C_SW,HCww.C_Eff,proj,HCww.OpTime,Cms,T,mu,\
                   0,Pt['MF'],Pt['MF'],C_Vant,1,inf,v,R,M_geometry,MW_NaCl)

roc_sw = pse.combo(membs,RPE,HCww_TBdpi,TBdp.Q_ROC,TBdp.C_ROC,TBdp.C_SW,proj,TBdp.OpTime,Cms,T,mu,\
                   0,Pt['MF'],Pt['MF'],C_Vant,2,inf,v,R,M_geometry,MW_NaCl)

#ROC-WW trial with Microfiltration enabled for WW side
#roc_ww3 = pse.combo(membs,RPE,SCRA_TBdpi,SCRA.Q_Eff*1.4,TBdp.C_ROC,SCRA.C_Eff,proj,TBdp.OpTime,Cms,T,mu,\
#                  1,Pt['MF'],Pt['MF'],C_Vant,2,inf,v,R,M_geometry,MW_NaCl)

#ROC-WW trial with Ultrafiltration enabled for WW side
#roc_ww4 = pse.combo(membs,RPE,SCRA_TBdpi,SCRA.Q_Eff*1.4,TBdp.C_ROC,SCRA.C_Eff,proj,TBdp.OpTime,Cms,T,mu,\
#                   1,Pt['UF'],Pt['UF'],C_Vant,2,inf,v,R,M_geometry,MW_NaCl)

#ROC-WW trial with hydroturbine
#roc_ww5 = pse.combo(membs,Hydro,SCRA_TBdpi,SCRA.Q_Eff*1.4,TBdp.C_ROC,SCRA.C_Eff,proj,TBdp.OpTime,Cms,T,mu,\
#                   1,Pt['UF'],Pt['UF'],C_Vant,0,inf,v,R,M_geometry,MW_NaCl)

#ROC-WW trial with no transport 
#roc_ww6 = pse.combo(membs,RPE,SCRA_TBdpi,SCRA.Q_Eff*1.4,TBdp.C_ROC,SCRA.C_Eff,proj,TBdp.OpTime,12.23,T,mu,\
#                   0,Pt['UF'],Pt['UF'],C_Vant,0,inf,v,R,M_geometry,MW_NaCl)

#mincost - based on unit cost
#comp_mincost = pse.Compcg(sw_ww,roc_ww,roc_sw,'Cost_Unit($/kW)','min')
comp_mincost2 = pse.Compcg(sw_ww,roc_ww2,roc_sw,'Cost_Unit($/kW)','min')

#maxWnet
#comp_maxWnet = pse.Compcg(sw_ww,roc_ww,roc_sw,'MW_net(MW)','max')




#%%#Plotting###################################################################
#different plot inputs below - needs revising, now proof of concept
#compplot: inputs are data, graph title
#compmc = plotter.CompPlot(comp_mincost,"Min Cost")
#compmwn = plotter.CompPlot(comp_maxWnet,"Max MW_net")
#close("all") #closes graphs from previous run

#StackPlot: imputs are: dataset,proj,title
#compstmc = plotter.StackPlotSE(comp_mincost2,"Min Cost - Pelton")
compstmc = plotter.StackPlotSE(comp_mincost2,"Min Cost - RPE")

#compstmw = plotter.StackPlotSE(comp_maxWnet,"Max Net Power")

compstmc = plotter.StackPlotCostUnit(comp_mincost2,proj,"Min Cost")
#compstmw = plotter.StackPlotCost(comp_maxWnet,"Max Net Power")

compstmc = plotter.StackPlotROIUnit(comp_mincost2,"Min Cost")
#compstmw = plotter.StackPlotROI(comp_maxWnet,"Max Net Power")


#ROITPlot: inputs are proj,comp (the comparison dataset),case (e.g. sw_ww),OpTime
ROIT_roc_ww = plotter.ROITPlot(proj,comp_mincost2,"roc_ww",TBdp.OpTime)

end = timer()
print(end - start) # Time in seconds

#ScatmPlot: inputs are: dataset1,dataset2,dataset3,x,y,averager,spec,spec_inc, line, SysCon2_inc,interact,title
#Note - interactive plot has title, non-interactive does not
df4 = plotter.ScatmPlot(sw_ww,roc_ww2,roc_sw,'SE_net(kWh/m^3)','PV_net($)','Q_de(L/hr)',comp_mincost2,1,0,0,0,'Net Specific Energy vs Unit Net Present Value for each tested membrane')
#df5 = plotter.ScatmPlot(sw_ww,roc_ww2,roc_sw,'SE_net(kWh/m^3)','PV_net($)','Q_de(L/hr)',comp_mincost2,1,0,0,1,'Net Specific Energy vs Unit Net Present Value for each tested membrane')
#compcpw = plotter.ScatmPlot(roc_sw,sw_ww,roc_ww,'SE_net(kWh/m^3)','n_protr','# of Process Trains vs Specific Energy Comparison for Tampa Bay')
#Note, if the JSON -serializable error comes up for mpld3, see this fix by ceprio: https://github.com/mpld3/mpld3/issues/441
