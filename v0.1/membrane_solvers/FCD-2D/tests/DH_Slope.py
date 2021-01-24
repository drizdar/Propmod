from formulas import RelativeDiffusivity, Alpha_Phi
import math
import json
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
T = 298
D = RelativeDiffusivity(T,0)
A_phi = Alpha_Phi(D, T)
assert A_phi == 0.3919478269650434, 'Alpha_Phi does not match value expected'
replaceNans = lambda n: np.nan if n is None else n
f = open('alpha_phi.json','r')
file = json.loads(f.read())
f.close()
actual = np.array(file['data'])
for i in range(0,len(actual)):
    for j in range(0,len(actual[i])):
        actual[i][j] = replaceNans(actual[i][j])
actual = np.array(actual, dtype=np.float32)
actual = np.transpose(actual)
P_arr = np.array(file['pressure'])
T_arrC = np.array(file['temperature'])
T_arrK = T_arrC + 273.15
predicted = np.zeros((len(P_arr),len(T_arrK)))
for i in range(0,len(P_arr)):
    for j in range(0,len(T_arrK)):
        P = P_arr[i]
        T = T_arrK[j]
        D = RelativeDiffusivity(T,P)
        predicted[i][j] = Alpha_Phi(D,T)

E1 = np.zeros((len(P_arr),len(T_arrK)))
E2 = np.zeros((len(P_arr),len(T_arrK)))
mean = np.nanmean(actual)
# for i in range(0,len(P_arr)):
for i in range(0,1):
    for j in range(0,len(T_arrK)):
        if (np.isnan(actual[i][j])):
            E1[i][j] = 0
            E2[i][j] = 0
        else:
            E1[i][j] = math.pow(actual[i][j] - predicted[i][j],2)
            E2[i][j] = math.pow(actual[i][j] - mean,2)
sumE1 = np.sum(E1)
sumE2 = np.sum(E2)
R2 = 1-sumE1/sumE2
print(R2)

fig = plt.figure(figsize=(9,5))
plt.axes([0.1,0.1,0.8,0.8], xlabel='Pressure - P (bar)', ylabel='Osmotic Coefficient - $Alpha_\phi$ (bar)')
plt.plot(T_arrK, actual[0,:], 'gs', T_arrK, predicted[0,:], 'r^')
plt.legend(['Observed', f'Model $R^2 = {R2:.5g}$', 'best'])
plt.show()