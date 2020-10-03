import sys
import math
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import formulas as f

###Test data
A = 0.67
B = 0.4
D = 6.25e-06
k = 306
S = 0.0008
n = 2
T = 25 + 273.15
R = f.R
Qdi = 18.181818181818183
Cdi = 3.6
Qfi = 47.61904761904762
Cfi = 0.6
dLd = 0.09236363636363637 
dLf = 0.09987042758895644
dA = dLd * dLf
H = 28*1e-3*0.0254
PId = f.OsP(Cdi, n, R, T)
PIf = f.OsP(Cfi, n, R, T)
dP = (PId - PIf)/2
Pdi = 2.8 + dP
Pfi = 2.8

def calculateElement(H, dLd, dLf, A, B, D, k, S, n, R, T, Cdi, Pdi, Qdi, Cfi, Pfi, Qfi):
    #Estimate Jw from inital data
    PId = f.OsP(Cdi, n, R, T)
    PIf = f.OsP(Cfi, n, R, T)
    dP = Pdi - Pfi
    dA = dLd * dLf
    try:
        Jw_est = f.WF(A, B, D, k, S, PId, PIf, dP, 10)
    except: 
        print("Jw calculation failed!")
        sys.exit(1)
    #Calculate other values
    Js = f.SF(B, D, k, S, Cdi, Cfi, Jw_est)
    Qdf = f.IDF(Qdi, Jw_est, dA)
    Cdf = f.IDC(Qdi, Cdi, Qdf, Js, dA)
    Qff = f.IFF(Qfi, Jw_est, dA)
    Cff = f.IFC(Qfi, Cfi, Qff, Js, dA)
    Qdb = (Qdi + Qdf)/2
    Cdb = (Cdi + Cdf)/2
    Qfb = (Qfi + Qff)/2
    Cfb = (Cfi + Cff)/2
    PId = f.OsP(Cdb, n, R, T)
    PIf = f.OsP(Cfb, n, R, T)

    #estimate Jw from bulk data
    try:
        Jw_post = f.WF(A, B, D, k, S, PId, PIf, dP, Jw_est)
    except: 
        print("Jw calculation failed!")
        sys.exit(1)

    #Converge results

    cr = Jw_est/Jw_post -1
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
        except: 
            print("Jw calculation failed!")
            sys.exit(1)
        cr = abs(Jw_pre/Jw_post -1)

    #Calculate final results
    Vdi = f.V(Qdi, dLf, H)
    Vdf = f.V(Qdf, dLf, H)
    Pdf = f.BPD(Pdi, Vdi, Vdf, 1400, 1400)
    Vfi = f.V(Qfi, dLd, H)
    Vff = f.V(Qff, dLd, H)
    Pff = f.BPD(Pfi, Vfi, Vff, 1200, 1200)
    #TODO other pressure losses like spacers and membrane friction
    #Print report to console
    # try: 
    assert (Jw_post == 7.7356835901890575), "Jw_post"
    assert (Js == 40.69055450754387), "Js"
    assert (Qdf == 18.253175189480984), "Qdf"
    assert (Cdf == 3.585574691788887), "Cdf"
    assert (Qff == 47.547690611384816), "Qdf" 
    assert (Cff == 0.6010355214799753), "Cff"
    assert (Pdf == 77.16868540275011), "Pdf"
    assert (Pff == 2.8000007285888513), "Pff"
    
    # except:
    #     print('test fail')
    return {"Jw": Jw_post, "Js": Js, "Cdf": Cdf, "Pdf": Pdf, "Qdf": Qdf, "Cff": Cff, "Pff": Pff, "Qff": Qff}
ret = calculateElement(H, dLd, dLf, A, B, D, k, S, n, R, T, Cdi, Pdi, Qdi, Cfi, Pfi, Qfi)
print(ret)
