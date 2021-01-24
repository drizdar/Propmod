import classes as cl
import formulas as f
import math
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys

def calculateElement(membrane, initial_draw, initial_feed, Jw_i):
    #Estimate Jw from inital data
    [d_h] = membrane.GetDimensions(["d_h"])
    [dA] = membrane.GetDimensions(["dAm"])
    [A, B, S] = membrane.GetProperties(["A", "B", "S"])
    D = f.D
    rho = initial_draw.data.get("density")
    pc_wt = initial_draw.data.get("pc_wt")
    vel_d = f.Velocity(membrane, initial_draw, "draw")
    vel_f = f.Velocity(membrane, initial_feed, "feed")
    k, k_LH = f.MassTransportCoefficient(D/60/60, d_h, pc_wt, rho*1000, vel_d)
    PId = initial_draw.data.get("PI")
    PIf = initial_feed.data.get("PI")
    Pd = initial_draw.data.get("P")
    Pf = initial_feed.data.get("P")
    dP = Pd-Pf
    try:
        Jw_est = f.WF(A, B, D, k_LH, S, PId, PIf, dP, Jw_i)
    except: 
        raise Exception("Jw calculation failed!")

    #Calculate other values
    Cdi = initial_draw.data.get("molar_concentration")
    Cfi = initial_feed.data.get("molar_concentration")
    Js = f.SF(B, D, k_LH, S, Cdi, Cfi, Jw_est)
    vel_f = f.Velocity(membrane, initial_feed, "feed")
    final_draw = f.IterateFlow(initial_draw, membrane, Js, Jw_est, "draw", vel_d)
    final_feed = f.IterateFlow(initial_feed, membrane, Js, Jw_est, "feed", vel_f)
    bulk_draw = cl.averageFlows(initial_draw, final_draw)
    bulk_feed = cl.averageFlows(initial_feed, final_feed)
    PId = initial_draw.data.get("PI")
    PIf = initial_feed.data.get("PI")
    rho = bulk_draw.data.get("density")
    pc_wt = bulk_draw.data.get("pc_wt")
    vel_d = f.Velocity(membrane, initial_draw, "draw")
    k, k_LH = f.MassTransportCoefficient(D/60/60, d_h, pc_wt, rho*1000, vel_d)
    Pd = bulk_draw.data.get("P")
    Pf = bulk_feed.data.get("P")
    Cd = bulk_draw.data.get("molar_concentration")
    Cf = bulk_feed.data.get("molar_concentration")
    dP = Pd-Pf
    #estimate Jw from bulk data

    try:
        Jw_post = f.WF(A, B, D, k_LH, S, PId, PIf, dP, Jw_est)
    except: 
        raise Exception("Jw calculation failed!")

    #Converge results

    cr = Jw_est/Jw_post -1
    while cr > 1e-12:
        Jw_pre = Jw_post
        Cdi = initial_draw.data.get("molar_concentration")
        Cfi = initial_feed.data.get("molar_concentration")
        Js = f.SF(B, D, k_LH, S, Cdi, Cfi, Jw_est)
        vel_f = f.Velocity(membrane, initial_feed, "feed")
        final_draw = f.IterateFlow(initial_draw, membrane, Js, Jw_pre, "draw", vel_d)
        final_feed = f.IterateFlow(initial_feed, membrane, Js, Jw_pre, "feed", vel_f)
        bulk_draw = cl.averageFlows(initial_draw, final_draw)
        bulk_feed = cl.averageFlows(initial_feed, final_feed)
        PId = initial_draw.data.get("PI")
        PIf = initial_feed.data.get("PI")
        rho = bulk_draw.data.get("density")
        pc_wt = bulk_draw.data.get("pc_wt")
        vel_d = f.Velocity(membrane, initial_draw, "draw")
        k, k_LH = f.MassTransportCoefficient(D/60/60, d_h, pc_wt, rho*1000, vel_d)
        Pd = bulk_draw.data.get("P")
        Pf = bulk_feed.data.get("P")
        Cd = bulk_draw.data.get("molar_concentration")
        Cf = bulk_feed.data.get("molar_concentration")
        dP = Pd-Pf

        try:
            Jw_post = f.WF(A, B, D, k_LH, S, PId, PIf, dP, Jw_pre)
        except: 
            raise Exception("Jw calculation failed!")
        cr = abs(Jw_pre/Jw_post -1)

    Pd = final_draw.data.get("P")
    return {"Cd": Cd, "Cf": Cf, "Jw": Jw_post, "Js": Js, "Pd": Pd, "Pf": Pf, "final_draw": final_draw, "final_feed": final_feed}