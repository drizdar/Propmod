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

def FeedBoundary(A, initial_feed, ned, nef, n_leafs):
    for j in range(0, ned): #rows, along feed direction
        A[0][j] = cl.splitFlows(initial_feed, ned*n_leafs)
    return A

def calculate(initial_draw, initial_feed):
    #assume one leaf per membrane
    n_PV = 2500
    n_leafs = 5
    Pi = 10 #bar
    initial_draw = cl.splitFlows(initial_draw, n_PV)
    initial_feed = cl.splitFlows(initial_feed, n_PV)
    fs = open('membrane_data.json', 'r')
    data = json.loads(fs.read())
    s = 0.1
    membrane = cl.membrane({
        "properties": data[0],
        "dimensions": {
            "Am": 40.08, #m^3
            "Ld": 1.016, #m
            "Hc": 0.0007, #m
            "Ss": 0.0045 #m
        }
    },s)
    [ned, nef] = membrane.GetDimensions(["ned", "nef"])
    Jw = np.zeros((nef, ned))
    Js = np.zeros((nef, ned))
    PPd = np.zeros((nef, ned))
    PId = initial_draw.data.get("PI")
    PIf = initial_feed.data.get("PI")
    Pd = initial_draw.data.get("P")
    Pf = initial_feed.data.get("P")
    dP = (PId-PIf)/2
    initial_draw.data["P"] += dP + Pi
    initial_feed.data["P"] += Pi
    draw = EmptyArray(ned, nef)
    feed = EmptyArray(ned, nef)
    draw = DrawBoundary(draw, initial_draw, ned, nef)
    feed = FeedBoundary(feed, initial_feed, ned, nef, n_leafs)
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
                PPd[i][j] = ret['Pd']
                draw[i][j+1] = ret['final_draw']
                feed[i+1][j] = ret['final_feed']
                Jw_pre = ret['Jw']

    Jw_avg = np.average(Jw)
    Js_avg = np.average(Js)
    P = np.average(PPd[:,-1])

    PDens = P*1e5*Jw_avg/60/60/1000

    return PDens

def graphJw(Jw, membrane):
    # levels = MaxNLocator(nbins=50).tick_values(Jw.min(), Jw.max())
    # instance which takes data values and translates those into levels.
    # cmap = plt.get_cmap('rainbow')
    # norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

    # fig, (ax0, ax1) = plt.subplots(nrows=2)
    # im = ax0.pcolormesh(x, y, Jw, cmap=cmap, norm=norm)
    # fig.colorbar(im, ax=ax0)
    # ax0.set_title('pcolormesh with levels')


    # # contours are *point* based plots, so convert our bound into point
    # # centers
    # # cf = ax1.contourf(x[:-1] + dLd/2., y[:-1] + dLf/2., Jw, levels=levels, cmap=cmap)

    # fig = plt.figure(figsize=(4,6),tight_layout=True)
    # ax = fig.add_axes([0.1,0.1,0.8,0.8])
    # cf = ax.contour(x, y, Jw, colors='black')
    # ax.clabel(cf, inline=True, fontsize=6)
    [dLd, dLf, Ld, Lf, ned, nef] = membrane.GetDimensions(["dLd", "dLf", "Ld", "Lf", "ned", "nef"])
    x = np.linspace(dLd/2, Ld-dLd/2, ned)
    y = np.linspace(Lf-dLf/2, dLf/2, nef)
    plt.figure()
    cp = plt.contour(x[:], y[:], Jw, colors='black')
    plt.clabel(cp, inline=True, fontsize=6)
    plt.title(r'Flux on membrane surface $(L\ h^{-1}\ m^{-2})$')
    plt.xlabel(r'Draw position $(m)$')
    plt.ylabel(r'Feed position $(m)$')
    plt.show()
    # plt.savefig('test.png')