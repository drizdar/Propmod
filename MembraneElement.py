import sys
import math
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import formulas as f

def calculateElement(H, dLd, dLf, A, B, D, k, S, n, R, T, Cdi, Pdi, Qdi, Cfi, Pfi, Qfi, Jw_i):
    #Estimate Jw from inital data
    PId = f.OsP(Cdi, n, R, T)
    PIf = f.OsP(Cfi, n, R, T)
    dP = Pdi - Pfi
    dA = dLd * dLf
    try:
        Jw_est = f.WF(A, B, D, k, S, PId, PIf, dP, Jw_i)
    except: 
        print("Jw calculation failed!")
        return None

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
        return None

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
            return None
        cr = abs(Jw_pre/Jw_post -1)

    #Calculate final results

    Vdi = f.V(Qdi, dLf, H)
    Vdf = f.V(Qdf, dLf, H)
    Pdf = f.BPD(Pdi, Vdi, Vdf, 1400, 1400)
    Vfi = f.V(Qfi, dLd, H)
    Vff = f.V(Qff, dLd, H)
    Pff = f.BPD(Pfi, Vfi, Vff, 1200, 1200)

    return {"Jw": Jw_post, "Js": Js, "Cdf": Cdf, "Pdf": Pdf, "Qdf": Qdf, "Cff": Cff, "Pff": Pff, "Qff": Qff}