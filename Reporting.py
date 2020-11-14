import matplotlib
import matplotlib.pyplot as plt
from matplotlib.dates import (MONTHLY, DateFormatter, rrulewrapper, RRuleLocator, drange)
import numpy as np
import pandas as pd
from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()


def reportInit():
    pd.set_option('display.max_rows', None)
    return np.array([["Value", "Description"],["-----", "-----------"]])
def pushReport(val, desc, rep):
    return np.append(rep, [[val, desc]], axis = 0)


def PlotTimeSeries(d, d_start, pond, list, x_label, y_labels, colours, legends):
    data = []
    for i in range(0, len(list)):
        item = list[i]
        data.append([])
        for p in pond:
            data[i].append(p.data.get(item))
    
    t = np.linspace(0, d-d_start, d-d_start+1)

    #Pond time-series
    axes = []
    fig, ax = plt.subplots()
    axes.append(ax)
    axes[0].set_xlabel(x_label)
    axes[0].set_ylabel(y_labels[0])
    lines = axes[0].plot(t, data[0], colours[0])
    for i in range(1, len(list)):
        axes.append(axes[-1].twinx())
        axes[-1].set_ylabel(y_labels[i])
        lines += axes[-1].plot(t, data[i], colours[i])
    
    plt.legend(lines, legends,loc=9) #or bbox_to_anchor=(0.2, 1.0)
    plt.tight_layout()
    plt.show()