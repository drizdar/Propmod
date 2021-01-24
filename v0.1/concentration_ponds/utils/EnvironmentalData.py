import datetime
import json
import numpy as np
import pandas as pd
from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.dates import (MONTHLY, DateFormatter, rrulewrapper, RRuleLocator, drange)

pan_evaporation_data = np.array([450, 350, 350, 200, 125, 80, 100, 150, 200, 300, 400, 400])
rainfall_data = np.array([17.5, 26.1, 9.8, 9.0, 7.0, 10.7, 4.3, 7.3, 12.2, 9.1, 13.5, 14.5])
def SmoothMonthlyData(data, rolling_avg):
    days_in_months = np.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])
    expanded_data = []
    for m in range(11, 12):
        month_value = data[m]
        days = days_in_months[m]
        for d in range(0, days):
            expanded_data = np.append(expanded_data, [month_value/days])
        for m in range(0, 12):
                month_value = data[m]
                days = days_in_months[m]
                for d in range(0, days):
                    expanded_data = np.append(expanded_data, [month_value/days])
        for m in range(0, 1):
                month_value = data[m]
                days = days_in_months[m]
                for d in range(0, days):
                    expanded_data = np.append(expanded_data, [month_value/days])
    transformed = []
    for d in range(31, 31+365):
        day_value = 0
        for n in range(d-rolling_avg,d+rolling_avg):
            day_value += expanded_data[n]/(rolling_avg*2)
        transformed = np.append(transformed, [day_value])
    return transformed

pan_evaporation = SmoothMonthlyData(pan_evaporation_data, rolling_avg=15)
rainfall = SmoothMonthlyData(rainfall_data, rolling_avg=15)
evap_factor = 0.7 #assume it's harder to evaporate concentrated salt solutions
evaporation = np.array([])
for i in pan_evaporation:
    evaporation = np.append(evaporation, [i*evap_factor])
permeability = 1e-9 #m/s
infiltration_loss = permeability*60*60*24/1000 #mm/day
infiltration_losses = []
for i in range(0, 365):
    infiltration_losses = np.append(infiltration_losses, [infiltration_loss])

fs = open('environmental_data.json', 'w')
fs.write(json.dumps({
    "rainfall": rainfall.tolist(), "evaporation": evaporation.tolist(), "infiltration_losses": infiltration_losses.tolist(),
    "rainfall_avg": np.average(rainfall), "evaporation_avg": np.average(evaporation), "infiltration_losses_avg": np.average(infiltration_losses)
}))

print(f'Average evaporation {np.average(evaporation)} mm/d')
print(f'Average infiltration losses {np.average(infiltration_losses)} mm/d')
print(f'Average rainfall {np.average(rainfall)} mm/d')

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
plt.savefig('./../capstone/images/env-gain-loss.png')
plt.show()