# -*- coding: utf-8 -*-
"""
Created on Tue Feb 5 21:23:17 2019

@author: Drizdar

Desalination energy recovery system designer Sensitivity and Uncertainty Analyzer

Note: Currently configured for SW-WW

Can choose between (in order of decreasing rigor):
    Saltelli Method
    Fast Fourier Transform Method
    Random Balanced Design - Fourier Amplitude Sensitivity Test
    
Note: For Y analysis, default is Cost_Unit($/kWh). For other Y variables, look at pro_sys_eqns.comboUSA for full list
"""

import pro_sys_eqns as pse
import TabularData as td
import pandas as pd
from SALib.sample import saltelli, latin, fast_sampler
from SALib.analyze import fast, rbd_fast, sobol
from timeit import default_timer as timer

#%%#Universal Input Parameters#################################################
csv_title = input('What will the test files be called? ')
csv_title_output = "".join((csv_title,'_output'))
csv_title_input = "".join((csv_title,'_input'))
csv_title_SI = "".join((csv_title,'_SI'))

method = input('Which Method will you use? Choices are FAST, RBD, or Sobol: ')
N_samples = int(input('How many samples should be generated per variable (Minimum value = 65)? '))
Notes = input('Write down any notes about this test: ')


start = timer() #used for timing the function

#C_Vant = 0.7345 #bar kg g-1 -> Van't Hoff Coefficient, from He, W., Wang, Y., Shaheed, M., 2015. Stand-alone seawater RO (reverse osmosis) desalination powered by PV (photovoltaic) and PRO (pressure retarded osmosis). Energy 86, 423–435.
v = 2 #Van't Hoff factor for NaCl, from "Lin, S., Straub, A., Elimelech, M., 2014. Thermodynamic limits of extractable energy by pressure retarded osmosis. Energy Environ Sci 7, 2706–2714."
R = 0.083144598 #L-bar/mol-K
MW_NaCl = 58.442769 #g/mol
C_Vant = (v*R/MW_NaCl) #L-bar/g-K - Van't Hoff Coefficient (Modified)
R = 2.30957E-06 #kWh/Mol-K #Change Units of R
g = 9.81 #m/s^2
inf = 1.2089 #inflation from September 2007 to November 2018 for cost estimates, from https://www.bls.gov/data/inflation_calculator.htm

#Note: Inflation value for Turbine (Found in pse) may also need to be adjusted depending on when model is ran
A_m = 40.88 #m^2 ~440 sq. ft, from VI email chain
L_m = 1.016 #m - membrane element length (A)
Mpv = 7 #number of membrane elements per pressure vessel
h_c = 0.0007 #m - channel height
l_f = 0.0045 #m - length between spacers


#%%Other Inputs################################################################
#Pretreatment data, spacer data, and N
Pt = td.Pt_dat

Pt_dat = pd.DataFrame.from_records([Pt[s].to_dict() for s, value in Pt.items()]) #DataFrame with all Pretreatment data - not in order
Pt_dat.set_index("ID", inplace=True)

sd = td.Spacer_dat
Spacers = pd.DataFrame.from_records([sd[s].to_dict(h_c,l_f) for s, value in sd.items()]) #values in to_dict are h_c, l_f, in m
Spacers.set_index("Shape", inplace=True)

shape = "1:3 ellipse"
M_geometry = [Spacers.loc[shape,"alpha"],Spacers.loc[shape,"beta"],1,Spacers.loc[shape,"d_h"],
              h_c, l_f, A_m, L_m, Mpv] #values are alpha, beta, gamma, d_h, h_c, l_f, A_m, L_m, Mpv

N_Yrs = 30 #Total # of years plant will run

#%%#Set Up Analysis############################################################

problem = {
    'num_vars': 19,
    'names': ['k', 'D', 'S', 'B', 'A', 'Turbeff', 'z1', 'z2', 'L','d', 'Q', 'C_D', 'C_F', 'EP', 'dEP', 'r', 'OpTime', 'Cms', 'T'],
    'bounds': [[76.7, 138.6],                   #k [L m-2 h-1] - Mass Transfer Coefficient
               [5.210e-6, 5.364e-6],            #D [m^2/h] - DIffusion Coefficient
               [1.49e-4, 9.87e-4],              #S [m] - Structural Parameter
               [0.08, 0.88],                    #B [L m-2 h-1] - Salt Permeability
               [1.23, 5.81],                    #A [L m-2 h-1 bar-1] - Water Permeability
               [0.85, 0.99],                    #Turbeff
               [2, 20],                         #z1 [m above sea level] - WW Plant Elevation
               [0, 4],                          #z2 [m above sea level] - RO Plant Elevation
               [0, 23500],                      #L [m]
               [0.1, 1],                        #d [m] 
               [0.2, 0.83],                     #Q [m^3/s] - Draw Flowrate (since V draw > V feed - usually)
               [60, 72],                        #C_D [g/kg] - Draw Concentration
               [27, 35],                        #C_F [g/kg] - Feed concentration
               [0.08, 0.34],                    #EP [$/kWh] - Energy price
               [0.003, 0.015],                  #dEP [$] - Change in energy price
               [0.01, 0.05],                    #r - Interest rate
               [0.6, 1],                        #OpTime - For RO Plant 
               [5.80, 14.18],                   #Cms [$/m^2]
               [14.5,31.6]]                     #T [deg C] - range from Tampa Bay Station 28 since Aug 2015
}

if method.lower() in ['fast']: 
    param_values = fast_sampler.sample(problem, N_samples) #For Fourier Amplitude Sensitivity Test
    print("FAST")
elif method.lower() in ['rbd']:
    param_values = latin.sample(problem, N_samples) #For RBD_FAST Method
    print("RBD")
elif method.lower() in ['sobol']:
    param_values = saltelli.sample(problem, N_samples) #For Sobol Method
    print("Sobol")

#Inputs are prob,N_Yrs,Tropp,PT_d,PT_f,C_Vant,PTopp,inf,v,R,A_m,L_m,Mpv,g,MW_NaCl,csv_name
#Tropp = 1 = transmission, 0 = no transmission
#PTopp: 0 - no pretreatment, 1 - Draw Pretreatment, 2 = Feed pretreatment, 3 = feed and draw pretreatment
Test = pse.comboUSA(param_values,N_Yrs,0,Pt['MF'],Pt['MF'],C_Vant,2,inf,v,R,M_geometry,g,MW_NaCl,csv_title_output)

#%%#Run Analysis###############################################################
Y = Test['PV_net($)'].values #Currently, Y is set for Net Present Value. Can be changed as needed

if method.lower() in ['fast']:
    Si = fast.analyze(problem, Y) #For Fourier Amplitude Sensitivity Test
elif method.lower() in ['rbd']:
    Si = rbd_fast.analyze(problem, Y, param_values) #For RBD-FAST Method
elif method.lower() in ['sobol']:
    Si = sobol.analyze(problem, Y) #For Sobol Method
    
    #Type: input the analysis type, e.g. "sobol"
#Si = pse.SA_indices(problem,Si,Y,method,csv_title_SI)
csv_title_SI_UC = "".join((csv_title_SI,'_PV_net')) #Configured For net present value - change as needed
Si_Var = pse.SA_indices(problem,Si,'PV_net($)',method,csv_title_SI_UC,Notes,N_samples)

#Placed at end so that files not unnecessarily generated for failed tests
sys_in = pse.SAin(problem,param_values,csv_title_input)

end = timer()
print(end - start) # Time in seconds
