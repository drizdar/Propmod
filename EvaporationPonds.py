import math
import pandas as pd
from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.dates import (MONTHLY, DateFormatter, rrulewrapper, RRuleLocator, drange)
import datetime
import formulas as f
import classes as cl
T = 298.15
P = 1.01325
pan_evaporation_data = [450, 350, 350, 200, 125, 80, 100, 150, 200, 300, 400, 400]
rainfall_data = [17.5, 26.1, 9.8, 9.0, 7.0, 10.7, 4.3, 7.3, 12.2, 9.1, 13.5, 14.5]
def SmoothMonthlyData(data, rolling_avg):
    days_in_months = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    expanded_data = []
    for m in range(11, 12):
        month_value = data[m]
        days = days_in_months[m]
        for d in range(0, days):
            expanded_data.append(month_value/days)
        for m in range(0, 12):
                month_value = data[m]
                days = days_in_months[m]
                for d in range(0, days):
                    expanded_data.append(month_value/days)
        for m in range(0, 1):
                month_value = data[m]
                days = days_in_months[m]
                for d in range(0, days):
                    expanded_data.append(month_value/days)
    transformed = []
    for d in range(31, 31+365):
        day_value = 0
        for n in range(d-rolling_avg,d+rolling_avg):
            day_value += expanded_data[n]/(rolling_avg*2)
        transformed.append(day_value)
    return transformed

pan_evaporation = SmoothMonthlyData(pan_evaporation_data, rolling_avg=15)
rainfall = SmoothMonthlyData(rainfall_data, rolling_avg=15)
evap_factor = 0.7 #assume it's harder to evaporate concentrated salt solutions
evaporation = []
for i in pan_evaporation:
    evaporation.append(i*evap_factor)

RO = {
    "recovery_rate": 0.4
}
RO_water = cl.solution({
    "pc_wt": 0.0001,
    "P": P,
    "T": T,
    "flow": 128 #MLD
})
seawater = cl.solution({
    "pc_wt": 0.035,
    "P": P,
    "T": T,
    "flow": RO_water.data["flow"] / RO["recovery_rate"] #MLD
})
seawater.CalcOsmoticProperties()
seawater.CalcMass()
RO_water.CalcOsmoticProperties()
RO_water.CalcMass()
brine_discharge = cl.solution({
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

A = 19.7*1e6 #m^2
D = 6000 #mm
m_brine_w = brine_discharge.data["mass_water"]
m_brine_NaCl = brine_discharge.data["mass_NaCl"]
rho_w = f.DensityWater(T)
def SizePond(A, D, evaporation, infiltration_loss, m_brine_w, m_brine_NaCl, rainfall, rho_w):
    h = 0
    while h < 6000:
        H = [cl.solution({
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
        while (h < 6000 and pc_wts[len(pc_wts)-1] < 0.265):
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
            soil_loss_water = A * infiltration_loss * rho * (1-pc_wt)
            seepage_NaCl = A * infiltration_loss * rho * pc_wt
            m_w += max(brine_water + rain - evap - soil_loss_water, 0)
            m_NaCl += max(brine_NaCl - soil_loss_water, 0)
            m_t = m_w + m_NaCl
            pond = cl.solution({
                "P": P,
                "T": T,
                "mass_water": m_w,
                "mass_NaCl": m_NaCl,
                "mass_total": m_t
            })
            pond.CalcPcWt()
            pc_wts = np.append(pc_wts, pond.data["pc_wt"])
            pond.CalcOsmoticProperties()
            pond.CalcFlow()
            H.append(pond)
            V = pond.data["flow"] / 1 #L
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
ax1.plot(t, hh, 'gs')
ax1.set_xlabel('Time - d (days)')
ax2 = ax1.twinx()
ax1.set_ylabel('Height of solution in pond - h (mm)')
ax2.plot(t, pc_wts, 'r')
ax2.set_ylabel('Concentration NaCl - (%wt)')
plt.show()

#Average losses and gains over the year
fig, ax = plt.subplots()
rule = rrulewrapper(MONTHLY, interval=1)
loc = RRuleLocator(rule)
formatter = DateFormatter('%b')
start_date = datetime.date(2021,1,1)
end_date = datetime.date(2022,1,1)
delta = datetime.timedelta(days=1)
dates = drange(start_date, end_date, delta)
plt.plot_date(dates, evaporation, 'r')
plt.plot_date(dates, rainfall, 'b')
plt.plot_date(dates, infiltration_losses, 'g')
ax.set_xlabel('Month')
ax.set_ylabel('Absolute Value - ($mm/m^2/day$)')
ax.xaxis.set_major_locator(loc)
ax.xaxis.set_major_formatter(formatter)
# ax.xaxis.set_tick_params(rotation=30)
plt.legend(['Evaporation', 'Rainfall', 'Soil Infiltration'])
plt.show()