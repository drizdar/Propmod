import numpy as np
import pandas as pd

def reportInit():
    pd.set_option('display.max_rows', None)
    return np.array([["Value", "Description"],["-----", "-----------"]])
def pushReport(val, desc, rep):
    return np.append(rep, [[val, desc]], axis = 0)