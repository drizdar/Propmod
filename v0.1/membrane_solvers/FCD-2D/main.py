import classes as cl
import EvaporationPonds as ep
import formulas as f
import json
import math
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.dates import (MONTHLY, DateFormatter, rrulewrapper, RRuleLocator, drange)
import MembraneFlow as mf
import numpy as np
import pandas as pd
from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()
import Reporting as r

T = 298.15
P = 1.01325

RO = {
    "recovery_rate": 0.4
}
RO_water = cl.flow({
    "pc_wt": 0.0001,
    "P": P,
    "T": T,
    "flow": 128e6 #L/d
})
seawater = cl.flow({
    "pc_wt": 0.035,
    "P": P,
    "T": T,
    "flow": RO_water.data["flow"] / RO["recovery_rate"] #MLD
})
seawater.CalcOsmoticProperties()
seawater.CalcMassRate()
RO_water.CalcOsmoticProperties()
RO_water.CalcMassRate()
RO_brine = cl.flow({
    "P": P,
    "T": T,
    "mass_water": seawater.data["mass_water"] - RO_water.data["mass_water"],
    "mass_NaCl": seawater.data["mass_NaCl"] - RO_water.data["mass_NaCl"],
    "mass_total": seawater.data["mass_water"] - RO_water.data["mass_water"] + \
        seawater.data["mass_NaCl"] - RO_water.data["mass_NaCl"]
})
RO_brine.CalcPcWt()
RO_brine.CalcOsmoticProperties()
RO_brine.CalcFlow()


D = 9000 #mm
A = 24.8e6 #m^2 play around with this until it looks good
d_start = 243 #start in September just as evaporation is increasing for Summer
d = d_start
n_years = 5

no_volume = {
    "A": A,
    "D": D,
    "P": P,
    "T": T,
    "mass_water": 0,
    "mass_NaCl": 0,
    "mass_total": 0,
    "level": 0,
    "pc_wt": 0,
    "mass_NaCl_solid": 450*A*2.17 #450mm bed existing
}
no_flow = {
    "P": P,
    "T": T,
    "mass_water": 0,
    "mass_NaCl": 0,
    "mass_total": 0,
    "pc_wt": 0,
    "flow": 0
}
pond = [cl.pond(no_volume)]
concentrate = [cl.flow(no_flow)]
discharge = [RO_brine]
no_PRO = {
    "J_w_avg": 0,
    "J_s_avg": 0,
    "Pd_avg": 0
}
PRO = [{
    "J_w_avg": 0,
    "J_s_avg": 0,
    "Pd_avg": 0
}]
mix_ratio = 1.209163317980615
Qb = RO_brine.GetFlow("MLD")
# This is where the logic for deciding whether to start PRO begins

while d < (d_start+n_years*365):
    pond = ep.IterateFlows(concentrate, d, D, discharge, pond)
    pond_pc_wt = pond[-1].data.get("pc_wt")
    pond_L = pond[-1].data.get("level")
    if (pond_L > 1000 and pond_pc_wt > 0.15):
        concentrate.append(cl.flow({
            "P": 1.01325,
            "T": 273.15 + 25,
            "pc_wt": pond_pc_wt,
            "flow": Qb*mix_ratio
        }))
        concentrate[-1].CalcOsmoticProperties()
        concentrate[-1].CalcMassRate()
        discharge.append(cl.combineFlows(discharge[0], concentrate[-1]))
    else:
        concentrate.append(cl.flow(no_flow))
        PRO.append(no_PRO)
        discharge.append(discharge[0])
    print(d, 
        pond[-1].data["level"],
        pond[-1].data["pc_wt"],
        concentrate[-1].data["flow"], 
        discharge[-1].data["flow"])
    d += 1

#Pond time-series
val_of_int = ["level", "pc_wt"]
y_labels = ['Depth of solution in pond - d (mm)','Concentration NaCl - (%wt)']
colours = ['g', 'r']
legends = ['Solution Depth', 'Concentration']
data = []
for i in range(0, len(val_of_int)):
    item = val_of_int[i]
    data.append([])
    for p in pond:
        data[i].append(p.data.get(item))
t = np.linspace(0, d-d_start, d-d_start+1)

axes = []
fig, ax = plt.subplots()
axes.append(ax)
axes[0].set_xlabel("Time - d (days)")
axes[0].set_ylabel(y_labels[0])
lines = axes[0].plot(t, data[0], colours[0])
for i in range(1, len(val_of_int)):
    axes.append(axes[-1].twinx())
    axes[-1].set_ylabel(y_labels[i])
    lines += axes[-1].plot(t, data[i], colours[i])

plt.legend(lines, legends,loc=9) #or bbox_to_anchor=(0.2, 1.0)
plt.tight_layout()
plt.savefig('./../capstone/images/pond-level-conc.png')
plt.show()





val_of_int = ["level", "level_NaCl"]
y_labels = ['Depth of solution in pond - d (mm)','Thickness of solid NaCl in pond - d (mm)']
colours = ['g', 'b']
legends = ['Solution Depth', 'Salt Layer Thickness']
data = []
for i in range(0, len(val_of_int)):
    item = val_of_int[i]
    data.append([])
    for p in pond:
        data[i].append(p.data.get(item))
t = np.linspace(0, d-d_start, d-d_start+1)

axes = []
fig, ax = plt.subplots()
axes.append(ax)
axes[0].set_xlabel("Time - d (days)")
axes[0].set_ylabel(y_labels[0])
lines = axes[0].plot(t, data[0], colours[0])
for i in range(1, len(val_of_int)):
    axes.append(axes[-1].twinx())
    axes[-1].set_ylabel(y_labels[i])
    lines += axes[-1].plot(t, data[i], colours[i])

plt.legend(lines, legends,loc=9) #or bbox_to_anchor=(0.2, 1.0)
plt.tight_layout()
plt.savefig('./../capstone/images/pond-salt-and-solution-level.png')
plt.show()