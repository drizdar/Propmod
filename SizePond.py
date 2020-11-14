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

Amax = 151*1e6 #m^2
T = 298.15
P = 1.01325
fs = open('environmental_data.json','r')
data = json.loads(fs.read())
rainfall = data.get("rainfall")
evaporation = data.get("evaporation")
infiltration_losses = data.get("infiltration_losses")

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
brine_discharge = cl.flow({
    "P": P,
    "T": T,
    "mass_water": seawater.data["mass_water"] - RO_water.data["mass_water"],
    "mass_NaCl": seawater.data["mass_NaCl"] - RO_water.data["mass_NaCl"],
    "mass_total": seawater.data["mass_water"] - RO_water.data["mass_water"] + \
        seawater.data["mass_NaCl"] - RO_water.data["mass_NaCl"]
})
brine_discharge.CalcPcWt()
brine_discharge.CalcOsmoticProperties()
brine_discharge.CalcFlow()

#Initial Stage
#Fill evaporation pond ready for PRO operation
permeability = 1e-9 #m/s
infiltration_loss = permeability*60*60*24/1000 #mm/day
infiltration_losses = []
for i in range(0, 365):
    infiltration_losses.append(infiltration_loss)

A = 21*1e6 #m^2
D = 9000 #mm
m_brine_w = brine_discharge.data["mass_water"]
m_brine_NaCl = brine_discharge.data["mass_NaCl"]
rho_w = f.DensityWater(T)

def SizePond(A, D, evaporation, infiltration_loss, m_brine_w, m_brine_NaCl, rainfall, rho_w):
    h = 0
    while h < D:
        H = [cl.flow({
                "P": P,
                "T": T,
                "mass_water": m_brine_w,
                "mass_NaCl": m_brine_NaCl,
                "mass_total": m_brine_w + m_brine_NaCl
            })]
        H[0].CalcPcWt()
        H[0].CalcOsmoticProperties()
        H[0].CalcFlow()
        d = 243 #1st September
        hh = np.array(H[0].data["flow"]/A)
        pc_wts = np.array([H[0].data["pc_wt"]])
        iterations = 0
        while (h < D and pc_wts[len(pc_wts)-1] < 0.265):
            d+=1
            previous_day = H[len(H)-1]
            dd = d % 365
            m_NaCl = previous_day.data["mass_NaCl"]
            m_w = previous_day.data["mass_water"]
            rho = previous_day.data.get("density")
            pc_wt = pc_wts[len(pc_wts)-1]
            brine_water = m_brine_w
            brine_NaCl = m_brine_NaCl
            rain = rainfall[dd] * A * rho_w #kg
            evap = evaporation[dd] * A * rho_w #kg
            soil_loss_water = A * infiltration_losses[dd] * rho * (1-pc_wt) if h > 0 else 0
            seepage_NaCl = A * infiltration_losses[dd] * rho * pc_wt if h > 0 else 0
            m_w += max(brine_water + rain - evap - soil_loss_water, 0)
            m_NaCl += max(brine_NaCl - soil_loss_water, 0)
            m_t = m_w + m_NaCl
            pond = cl.pond({
                "P": P,
                "T": T,
                "mass_water": m_w,
                "mass_NaCl": m_NaCl,
                "mass_total": m_t
            })
            pond.CalcPcWt()
            pc_wts = np.append(pc_wts, pond.data["pc_wt"])
            pond.CalcOsmoticProperties()
            pond.CalcVolume()
            H.append(pond)
            V = pond.data["V"] / 1 #L
            h = V/A
            hh = np.append(hh, h)
            iterations += 1
        A *= 0.95
        print(h, iterations, A)
        return [d, hh, iterations, pc_wts]
d, hh, iterations, pc_wts = SizePond(A, D, evaporation, infiltration_loss, m_brine_w, m_brine_NaCl, rainfall, rho_w)
t = np.linspace(0, d-243, d+1-243)

#Pond time-series
fig1, ax1 = plt.subplots()
ln1 = ax1.plot(t, hh, 'g')
ax1.set_xlabel('Time - d (days)')
ax2 = ax1.twinx()
ax1.set_ylabel('Height of solution in pond - h (mm)')
ln2 = ax2.plot(t, pc_wts, 'r')
ax2.set_ylabel('Concentration NaCl - (%wt)')
plt.legend(ln1+ln2,['Solution Level', 'Concentration'],loc=9) #or bbox_to_anchor=(0.2, 1.0)
plt.tight_layout()
# ax2.legend(['Concentration'], loc='upper right')
# plt.legend(['Solution Level', 'Concentration'])
plt.savefig('./../capstone/images/pond-timeseries.png')
plt.show()
