import classes as cl
import datetime
import formulas as f
import json
import math
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.dates import (MONTHLY, DateFormatter, rrulewrapper, RRuleLocator, drange)
import numpy as np
import pandas as pd
from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()

fs = open('environmental_data.json','r')
data = json.loads(fs.read())
rainfall = data.get("rainfall")
evaporation = data.get("evaporation")
infiltration_losses = data.get("infiltration_losses")



#Initial Stage
#Fill evaporation pond ready for PRO operation
def IterateFlows(brine_discharge, pond, d):
    d+=1
    dd = d % 365
    A = pond[-1].data.get("A")
    P = brine_discharge[-1].data.get("P")
    T = brine_discharge[-1].data.get("T")
    m_NaCl = pond[-1].data["mass_NaCl"]
    m_w = pond[-1].data["mass_water"]
    rho = pond[-1].data.get("density")
    rho_w = f.DensityWater(T)
    pc_wt = pond[-1].data.get("pc_wt")
    brine_water = brine_discharge[-1].data["mass_water"]
    brine_NaCl = brine_discharge[-1].data["mass_NaCl"]
    rain = rainfall[dd] * A * rho_w #kg
    evap = evaporation[dd] * A * rho_w #kg
    soil_loss_water = A * infiltration_losses[dd] * rho * (1-pc_wt) if pond[-1].data["level"] > 0 else 0
    seepage_NaCl = A * infiltration_losses[dd] * rho * pc_wt if pond[-1].data["level"] > 0 else 0
    m_w += max(brine_water + rain - evap - soil_loss_water, 0)
    m_NaCl += max(brine_NaCl - soil_loss_water, 0)
    m_t = m_w + m_NaCl
    pond.append(cl.pond({
        "A": A,
        "P": P,
        "T": T,
        "mass_water": m_w,
        "mass_NaCl": m_NaCl,
        "mass_total": m_t,
    }))
    pond[-1].CalcPcWt()
    pond[-1].CalcOsmoticProperties()
    pond[-1].CalcVolume()
    pond[-1].CalcLevel()
    return [brine_discharge, pond]