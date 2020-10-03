import sys
import math
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import formulas as f
def reportInit():
    pd.set_option('display.max_rows', None)
    return np.array([["Value", "Description"]])
def pushReport(val, desc, rep): 
    return np.append(rep, [[val, desc]], axis = 0)
reportingOn = True
rep = reportInit() if reportingOn else None
rep = pushReport("","----Test data", rep) if reportingOn else rep
###Test data
A = 0.67; rep = pushReport(A, "A: Water permeability coefficient (L/h/m^2/bar)", rep) if reportingOn else rep
B = 0.4; rep = pushReport(B, "B: Salt permeability coefficient (g/h/m^2)", rep) if reportingOn else rep
D = 6.25e-06; rep = pushReport(D, "D: Diffusion coefficient for NaCL (m^2/h)", rep) if reportingOn else rep
k = 306; rep = pushReport(k, "k: External concentration polarisation coefficient", rep) if reportingOn else rep #TODO ADD UNTIS
S = 0.0008; rep = pushReport(S, "S: Structural parameter", rep) if reportingOn else rep
T = 25 + 273.15; rep = pushReport(T, "T: Temperature (K)", rep) if reportingOn else rep
R = f.R; rep = pushReport(R, "R: Ideal Gas Constant (mol/bar/K)", rep) if reportingOn else rep #TODO DOUBLE CHECK UNITS
Qdi = 18.181818181818183; rep = pushReport(Qdi, "Qdi: Initial draw flow across element (L/h)", rep) if reportingOn else rep
Cdi = 3.6; rep = pushReport(Cdi, "Cdi: Initial draw concentration (mol/L)", rep) if reportingOn else rep
Qfi = 47.61904761904762; rep = pushReport(Qfi, "Qfi: Initial feed flow across element (L/h)", rep) if reportingOn else rep
Cfi = 0.6; rep = pushReport(Cfi, "Cfi: Initial feed concentration (mol/L)", rep) if reportingOn else rep
PId = f.OsP(Cdi, 2, R, T); rep = pushReport(PId, "PId: Draw Osmotic Pressure (bar)", rep) if reportingOn else rep
PIf = f.OsP(Cfi, 2, R, T); rep = pushReport(PIf, "PIf: Feed Osmotic Pressure (bar)", rep) if reportingOn else rep
dP = (PId - PIf)/2; rep = pushReport(dP, "dP: Inital optimum pressure (bar)", rep) if reportingOn else rep
dLd = 0.09236363636363637 
dLf = 0.09987042758895644
dA = dLd * dLf; rep = pushReport(dA, "dA: Area of element (m^2)", rep) if reportingOn else rep
H = 28*1e-3*0.0254; rep = pushReport(H,'H: Channel height (m)', rep) if reportingOn else rep #SW30XHR-400i assume same height for both
Pdi = 2.8 + dP; rep = pushReport(Pdi,'Pdi: Initial draw pressure (bar)', rep) if reportingOn else rep
Pfi = 2.8; rep = pushReport(Pfi,'Pfi: Intial feed pressure (bar)', rep) if reportingOn else rep

###Function space
#Estimate Jw from inital data
rep = pushReport("", "----Estimate Jw from intial data", rep) if reportingOn else rep
try:
    Jw_est = f.WF(A, B, D, k, S, PId, PIf, dP, 10); rep = pushReport(Jw_est, "Jw_est: Water flux (L/h/m^2/bar)", rep) if reportingOn else rep
    print(Jw_est)
except: 
    print("Jw calculation failed!")
    sys.exit(1)
#Calculate other values
Js = f.SF(B, D, k, S, Cdi, Cfi, Jw_est); rep = pushReport(Js, "Js: Salt flux (g/h/m^2)", rep) if reportingOn else rep
Qdf = f.IDF(Qdi, Jw_est, dA); rep = pushReport(Qdf, "Qdf: Final draw flow across element (L/h)", rep) if reportingOn else rep
Cdf = f.IDC(Qdi, Cdi, Qdf, Js, dA); rep = pushReport(Cdf, "Cdf: Final draw concentration (mol/L)", rep) if reportingOn else rep
Qff = f.IFF(Qfi, Jw_est, dA); rep = pushReport(Qff, "Qdf: Final feed flow across element (L/h)", rep) if reportingOn else rep
Cff = f.IFC(Qfi, Cfi, Qff, Js, dA); rep = pushReport(Cff, "Cff: Final feed concentration (mol/L)", rep) if reportingOn else rep
Qdb = (Qdi + Qdf)/2; rep = pushReport(Qdb, "Qdb: Bulk draw flow across membrane (L/h)", rep) if reportingOn else rep
Cdb = (Cdi + Cdf)/2; rep = pushReport(Cdb, "Cdb: Bulk draw concentration (mol/L)", rep) if reportingOn else rep
Qfb = (Qfi + Qff)/2; rep = pushReport(Qfb, "Qfb: Bulk feed flow across element (L/h)", rep) if reportingOn else rep
Cfb = (Cfi + Cff)/2; rep = pushReport(Cfb, "Cfb: Bulk feed concentration (mol/L)", rep) if reportingOn else rep
PId = f.OsP(Cdb, 2, R, T); rep = pushReport(PId,'PId: Draw Osmotic Pressure (bar)', rep) if reportingOn else rep
PIf = f.OsP(Cfb, 2, R, T); rep = pushReport(PIf,'PIf: Feed Osmotic Pressure (bar)', rep) if reportingOn else rep
rep = pushReport("", "----Estimate Jw from bulk results", rep) if reportingOn else rep
#estimate Jw from bulk data
try:
    Jw_post = f.WF(A, B, D, k, S, PId, PIf, dP, Jw_est); rep = pushReport(Jw_post, "Jw_post: Water flux (L/h/m^2/bar)", rep) if reportingOn else rep
    print(Jw_post)
except: 
    print("Jw calculation failed!")
    sys.exit(1)

#Converge results
rep = pushReport("","----Concerge results", rep) if reportingOn else rep
cr = Jw_est/Jw_post -1; rep = pushReport(cr,'cr: Covergence Ratio', rep) if reportingOn else rep
while cr > 1e-12:
    Jw_pre = Jw_post
    Js = f.SF(B, D, k, S, Cdi, Cfi, Jw_pre)
    Qdf = f.IDF(Qdi, Jw_pre, dA)
    Qff = f.IFF(Qfi, Jw_pre, dA)
    Cff = f.IFC(Qfi, Cfi, Qff, Js, dA)
    Qdb = (Qdi + Qdf)/2
    Cdb = (Cdi + Cdf)/2
    Qfb = (Qfi + Qff)/2
    Cfb = (Cfi + Cff)/2
    PId = f.OsP(Cdb, 2, R, T)
    PIf = f.OsP(Cfb, 2, R, T)
    Cdf = f.IDC(Qdi, Cdi, Qdf, Js, dA)
    try:
        Jw_post = f.WF(A, B, D, k, S, PId, PIf, dP, Jw_pre)
        print(Jw_post)
    except: 
        print("Jw calculation failed!")
        sys.exit(1)
    cr = abs(Jw_pre/Jw_post -1); rep = pushReport(Jw_pre/Jw_post -1,'cr: Covergence Ratio', rep) if reportingOn else rep

#Calculate final results
rep = pushReport("","----Converged final results", rep) if reportingOn else rep
rep = pushReport(Jw_post, "Jw_final: Water flux (L/h/m^2/bar)", rep) if reportingOn else rep
rep = pushReport(Js, "Js: Salt flux (g/h/m^2)", rep) if reportingOn else rep
rep = pushReport(Qdf, "Qdf: Final draw flow across element (L/h)", rep) if reportingOn else rep
rep = pushReport(Cdf, "Cdf: Final draw concentration (mol/L)", rep) if reportingOn else rep
rep = pushReport(Qff, "Qdf: Final feed flow across element (L/h)", rep) if reportingOn else rep
rep = pushReport(Cff, "Cff: Final feed concentration (mol/L)", rep) if reportingOn else rep
rep = pushReport(Qdb, "Qdb: Bulk draw flow across membrane (L/h)", rep) if reportingOn else rep
rep = pushReport(Cdb, "Cdb: Bulk draw concentration (mol/L)", rep) if reportingOn else rep
rep = pushReport(Qfb, "Qfb: Bulk feed flow across element (L/h)", rep) if reportingOn else rep
rep = pushReport(Cfb, "Cfb: Bulk feed concentration (mol/L)", rep) if reportingOn else rep
rep = pushReport(PId,'PId: Draw Osmotic Pressure (bar)', rep) if reportingOn else rep
rep = pushReport(PIf,'PIf: Feed Osmotic Pressure (bar)', rep) if reportingOn else rep
Vdi = f.V(Qdi, dLf, H); rep = pushReport(Vdi,'Vdi: Initial draw velocity (m/s)', rep) if reportingOn else rep #TODO drop this part, it's almost negligible
Vdf = f.V(Qdf, dLf, H); rep = pushReport(Vdf,'Vdf: Final draw velocity (m/s)', rep) if reportingOn else rep
Pdf = f.BPD(Pdi, Vdi, Vdf, 1400, 1400); rep = pushReport(Pdf,'Pdf: Final draw pressure (bar)', rep) if reportingOn else rep #TODO low priority calculate density change
Vfi = f.V(Qfi, dLd, H); rep = pushReport(Vfi,'Vdi: Initial feed velocity (m/s)', rep) if reportingOn else rep
Vff = f.V(Qff, dLd, H); rep = pushReport(Vff,'Vdf: Final feed velocity (m/s)', rep) if reportingOn else rep
Pff = f.BPD(Pfi, Vfi, Vff, 1200, 1200); rep = pushReport(Pff,'Pff: Final feed pressure (bar)', rep) if reportingOn else rep #TODO low priority calculate density change
#TODO other pressure losses like spacers and membrane friction
#Print report to console
if reportingOn: 
    print(pd.DataFrame(rep[:,0],rep[:,1],[""]))