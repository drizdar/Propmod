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

dissolution_rate = 2.5 #mm/d
fs = open('environmental_data.json','r')
data = json.loads(fs.read())
rainfall = data.get("rainfall")
evaporation = data.get("evaporation")
infiltration_losses = data.get("infiltration_losses")

def IterateFlows(concentrate, d, D, discharge, pond): 
    d+=1
    dd = d % 365
    if concentrate[-1].data.get("flow") > 0: 
        flow = concentrate[-1].data.get("flow")
        concentrate[-1].data["pc_wt"] = pond[-1].data.get("pc_wt")
        concentrate[-1].data["P"] = pond[-1].data.get("P")
        concentrate[-1].data["T"] = pond[-1].data.get("T")
        concentrate[-1].CalcOsmoticProperties()
        concentrate[-1].CalcMassRate()
        concentrate_water = concentrate[-1].data.get("mass_water")
        concentrate_NaCl = concentrate[-1].data.get("mass_NaCl")
    else: 
        concentrate_water = concentrate[-1].data.get("mass_water")
        concentrate_NaCl = concentrate[-1].data.get("mass_NaCl")
    A = pond[-1].data.get("A")
    P = discharge[-1].data.get("P")
    T = discharge[-1].data.get("T")
    m_NaCl = pond[-1].data["mass_NaCl"]
    m_w = pond[-1].data["mass_water"]
    rho = pond[-1].data.get("density")
    rho_w = f.DensityWater(T)
    pc_wt = pond[-1].data.get("pc_wt")
    discharge_water = discharge[-1].data["mass_water"]
    discharge_NaCl = discharge[-1].data["mass_NaCl"]
    rain = rainfall[dd] * A * rho_w #kg
    evap = evaporation[dd] * A * rho_w #kg
    soil_loss_water = A * infiltration_losses[dd] * rho * (1-pc_wt) if pond[-1].data["level"] > 0 else 0
    seepage_NaCl = A * infiltration_losses[dd] * rho * pc_wt if pond[-1].data["level"] > 0 else 0
    d_m_w = discharge_water + rain - evap - soil_loss_water - concentrate_water
    m_w += d_m_w
    if m_w < 0: raise Exception("All water lost")
    d_m_NaCl = discharge_NaCl - soil_loss_water - concentrate_NaCl
    m_NaCl += d_m_NaCl
    if m_NaCl < 0: raise Exception("All NaCl lost")
    m_t = m_w + m_NaCl
    m_NaCl_s = pond[-1].data.get("mass_NaCl_solid")
    pond.append(cl.pond({
        "A": A,
        "P": P,
        "T": T,
        "mass_water": m_w,
        "mass_NaCl": m_NaCl,
        "mass_total": m_t,
        "mass_NaCl_solid": m_NaCl_s,
    }))
    pond[-1].CalcPcWt()
    pond[-1].CalcNaClSolidLevel()
    pond_pc_wt = pond[-1].data.get("pc_wt")
    if pond_pc_wt > 0.265: #saturation level
        pond_m_t = pond[-1].data.get("mass_total")
        pond_m_NaCl = pond[-1].data.get("mass_NaCl")
        NaCl_to_solid = pond_m_t*(pond_pc_wt-0.265)
        pond[-1].data["mass_NaCl_solid"] += NaCl_to_solid
        pond[-1].data["mass_NaCl"] -= NaCl_to_solid
    else:
        m_NaCl_s = pond[-1].data.get("mass_NaCl_solid")
        if m_NaCl_s > 0:
            pond_m_t = pond[-1].data.get("mass_total")
            NaCl_to_aqueous = min(pond[-1].CalcMassOfDepthSolidNaCl(dissolution_rate), pond_m_t*(0.265-pond_pc_wt))
            pond[-1].data["mass_NaCl_solid"] -= NaCl_to_aqueous
            pond[-1].data["mass_NaCl"] += NaCl_to_aqueous
    pond[-1].CalcPcWt()
    pond[-1].CalcOsmoticProperties()
    pond[-1].CalcVolume()
    pond[-1].CalcLevel()
    pond[-1].CalcNaClSolidLevel()
    cr = 1
    post = pond[-1].data["mass_NaCl"]
    while cr > 1e-12:
        pre = post
        pond_pc_wt = pond[-1].data.get("pc_wt")
        if pond_pc_wt > 0.265: #saturation level
            pond_m_t = pond[-1].data.get("mass_total")
            pond_m_NaCl = pond[-1].data.get("mass_NaCl")
            NaCl_to_solid = pond_m_t*(pond_pc_wt-0.265)
            pond[-1].data["mass_NaCl_solid"] += NaCl_to_solid
            pond[-1].data["mass_NaCl"] -= NaCl_to_solid
        else:
            m_NaCl_s = pond[-1].data.get("mass_NaCl_solid")
            if m_NaCl_s > 0:
                pond_m_t = pond[-1].data.get("mass_total")
                NaCl_to_aqueous = min(pond[-1].CalcMassOfDepthSolidNaCl(dissolution_rate), pond_m_t*(0.265-pond_pc_wt))
                pond[-1].data["mass_NaCl_solid"] -= NaCl_to_aqueous
                pond[-1].data["mass_NaCl"] += NaCl_to_aqueous
        post = pond[-1].data["mass_NaCl"]
        pond[-1].CalcPcWt()
        pond[-1].CalcOsmoticProperties()
        pond[-1].CalcVolume()
        pond[-1].CalcLevel()
        pond[-1].CalcNaClSolidLevel()
        cr = abs(pre/post -1)

    return pond