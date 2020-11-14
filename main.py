import classes as cl
import EvaporationPonds as ep
import formulas as f
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
gamma_D = 0.3
gamma_F = 1-gamma_D
Qb = RO_brine.GetFlow("L/d")
Qc = Qb * gamma_D/gamma_F

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
    "mass_NaCl_solid": 0
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
            "flow": Qc
        }))
        concentrate[-1].CalcOsmoticProperties()
        concentrate[-1].CalcMassRate()
        discharge.append(cl.combineFlows(discharge[0], concentrate[-1]))
        PRO.append(mf.calculate(RO_brine, concentrate))
    else:
        concentrate.append(cl.flow(no_flow))
        PRO.append(no_PRO)
        discharge.append(discharge[0])
    print(d, 
        pond[-1].data["level"],
        pond[-1].data["pc_wt"],
        concentrate[-1].data["flow"], 
        discharge[-1].data["flow"],
        pond[-1].data["level_NaCl"],
        pond[-1].data["mass_NaCl_solid"])
    d += 1


r.PlotTimeSeries(d, d_start, pond, 
    list=["level", "pc_wt"],
    x_label="Time - d (days)",
    y_labels=['Depth of solution in pond - d (mm)','Concentration NaCl - (%wt)'],
    colours=['g', 'r'],
    legends=['Solution Depth', 'Concentration'])

r.PlotTimeSeries(d, d_start, pond, 
    list=["level", "level_NaCl"],
    x_label="Time - d (days)",
    y_labels=['Depth of solution in pond - d (mm)','Thickness of solid NaCl in pond - d (mm)'],
    colours=['g', 'b'],
    legends=['Solution Depth', 'Salt Layer Thickness'])

r.PlotTimeSeries(d, d_start, concentrate, 
    list=["flow", "pc_wt"],
    x_label="Time - d (days)",
    y_labels=['Flow (L/d)','Concentration NaCl - (%wt)'],
    colours=['g', 'b'],
    legends=['Pond Draw Solution', 'Concentration'])

r.PlotTimeSeries(d, d_start, discharge, 
    list=["flow", "pc_wt"],
    x_label="Time - d (days)",
    y_labels=['Flow (L/d)','Concentration NaCl - (%wt)'],
    colours=['g', 'b'],
    legends=['PRO Discharge', 'Concentration'])