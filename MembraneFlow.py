import classes as cl
import formulas as f
import json
import math
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import MembraneElement as me
import numpy as np
import pandas as pd
import Reporting as r
import sys

def EmptyArray(ned, nef):
    A = []
    for i in range(0, nef): #rows, along feed direction
        A.append([])
        for j in range(0,ned): #columns, along draw direction
            A[i].append([])
    return A

def DrawBoundary(A, initial_draw, ned, nef):
    for i in range(0, nef): #rows, along feed direction
        A[i][0] = cl.splitFlows(initial_draw, nef)
    return A

def FeedBoundary(A, initial_feed, ned, nef, n_leaves):
    for j in range(0, ned): #rows, along feed direction
        A[0][j] = cl.splitFlows(initial_feed, ned*n_leaves)
    return A

def calculate(membrane, initial_draw, initial_feed):
    [Am, ned, nef] = membrane.GetDimensions(["Am", "ned", "nef"])
    [n_leaves, n_pv] = membrane.GetProperties(["n_leaves", "n_pv"])
    Jw = np.zeros((nef, ned))
    Js = np.zeros((nef, ned))
    PPd = np.zeros((nef, ned))
    PPf = np.zeros((nef, ned))
    Cd = np.zeros((nef, ned))
    Cf = np.zeros((nef, ned))
    Qdi = initial_draw.GetFlow("L/d")/24 
    Qfi = initial_feed.GetFlow("L/d")/24 
    PId = initial_draw.data.get("PI")
    PIf = initial_feed.data.get("PI")
    Pi = 10 #intial pressure for both streams into PRO plant, bar
    Pd = initial_draw.data.get("P")
    Pf = initial_feed.data.get("P")
    dP = (PId-PIf)/2
    initial_draw.data["P"] += dP + Pi
    initial_feed.data["P"] += Pi
    pv_draw = cl.splitFlows(initial_draw, n_pv)
    pv_feed = cl.splitFlows(initial_feed, n_pv)
    draw = EmptyArray(ned+1, nef)
    feed = EmptyArray(ned, nef+1)
    draw = DrawBoundary(draw, pv_draw, ned+1, nef)
    feed = FeedBoundary(feed, pv_feed, ned, nef+1, n_leaves)
    Jw_pre = 10


    for i in range(0, nef): #rows, along feed direction
        for j in range(0,ned): #columns, along draw direction
            ret = me.calculateElement(membrane, draw[i][j], feed[i][j], Jw_pre)
            if ret == None:
                print('exiting')
                sys.exit(1)
            else:
                Jw[i][j] = ret['Jw']
                Js[i][j] = ret['Js']
                Cd[i][j] = ret['Cd']
                Cf[i][j] = ret['Cf']
                PPd[i][j] = ret['Pd']
                PPf[i][j] = ret['Pf']
                draw[i][j+1] = ret['final_draw']
                feed[i+1][j] = ret['final_feed']
                Jw_pre = ret['Jw']

    Jw_avg = np.average(Jw)
    Js_avg = np.average(Js)
    Pd_out = np.average(PPd[:,-1])
    Pf_out = np.average(PPf[-1,:])

    PDens = Pd_out*1e5*Jw_avg/60/60/1000
    DeltaP = Pd_out - Pf_out
    DeltaW = Jw_avg * Am
    SE = DeltaW * DeltaP / (Qdi+Qfi)

    return [Cd, Cf, Jw, Jw_avg, Js, Js_avg, PPd, PPf, PDens, SE]
