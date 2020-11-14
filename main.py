import classes as cl
import EvaporationPonds as ep
import formulas as f
import math
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.dates import (MONTHLY, DateFormatter, rrulewrapper, RRuleLocator, drange)
import numpy as np
import pandas as pd
from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()

T = 298.15
P = 1.01325

RO = {
    "recovery_rate": 0.4
}
RO_water = cl.flow({
    "pc_wt": 0.0001,
    "P": P,
    "T": T,
    "flow": 128 #MLD
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
brine_discharge = [cl.flow({
    "P": P,
    "T": T,
    "mass_water": seawater.data["mass_water"] - RO_water.data["mass_water"],
    "mass_NaCl": seawater.data["mass_NaCl"] - RO_water.data["mass_NaCl"],
    "mass_total": seawater.data["mass_water"] - RO_water.data["mass_water"] + \
        seawater.data["mass_NaCl"] - RO_water.data["mass_NaCl"]
})]
brine_discharge[-1].CalcPcWt()
brine_discharge[-1].CalcOsmoticProperties()
brine_discharge[-1].CalcFlow()

# This is where the logic for deciding whether to start PRO begins

D = 9000 #mm
A = 19.7*1e6 #m^2 from SizePond.py
pond = [cl.pond({
    "A": A,
    "D": D,
    "P": P,
    "T": T,
    "mass_water": 0,
    "mass_NaCl": 0,
    "mass_total": 0,
    "level": 0,
    "pc_wt": 0
})]
d = 243 #start in September just as evaporation is increasing for Summer

#Filling
while d < (243+4*365):
    brine_discharge, pond = ep.IterateFlows(brine_discharge, pond, d, D)
    print(pond[-1].data["level"], pond[-1].data["pc_wt"])
    d += 1
level = []
pc_wts = []
for p in pond:
    level.append(p.data["level"])
    pc_wts.append(p.data["pc_wt"])
t = np.linspace(0, d-243, d+1-243)
#Remove first record, it makes the graph look weird
level.pop()
pc_wts.pop()
t = t[:-1]
#Pond time-series
fig1, ax1 = plt.subplots()
ln1 = ax1.plot(t, level, 'g')
ax1.set_xlabel('Time - d (days)')
ax2 = ax1.twinx()
ax1.set_ylabel('Height of solution in pond - h (mm)')
ln2 = ax2.plot(t, pc_wts, 'r')
ax2.set_ylabel('Concentration NaCl - (%wt)')
plt.legend(ln1+ln2,['Solution Level', 'Concentration'],loc=9) #or bbox_to_anchor=(0.2, 1.0)
plt.tight_layout()
plt.show()

#Operating
