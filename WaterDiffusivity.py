from formulas import RelativeDiffusivity
import math
import json
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator

f = open('diffusivity.json','r')
data = json.loads(f.read())
actual = np.asarray(data)

T = [298]
P = actual[:,0]
predicted = np.zeros((len(T),len(P)))
for i in range(0,len(T)):
    for j in range(0,len(P)):
        predicted[i][j] = RelativeDiffusivity(T[i],P[j])

E1 = np.zeros(len(P))
E2 = np.zeros(len(P))
mean = np.mean(actual[:,1])
for i in range(0, len(predicted[0])):
    E1[i] = math.pow(actual[i][1] - predicted[0][i],2)
    E2[i] = math.pow(actual[i][1] - mean,2)
sumE1 = np.sum(E1)
sumE2 = np.sum(E2)
R2 = 1-sumE1/sumE2
print(R2)

assert R2 == 0.997460771783034, 'R2 value does not match expected'

fig = plt.figure(figsize=(9,5))
plt.axes([0.1,0.1,0.8,0.8], xlabel='Pressure - P (bar)', ylabel='Relative Diffusivity - D')
plt.plot(P, actual[:,1], 'gs', P, predicted[0], 'r^')
plt.legend(['Observed', f'Model $R^2 = {R2:.5g}$', 'best'])
plt.show()

# T = [298,338,373,523,573,623,673]
# P = np.linspace(0,4000,200)
# R = np.zeros((len(T),len(P)))
# for i in range(0,len(T)):
#     for j in range(0,len(P)):
#         C = U4 + U5/(U6 + T[i])
#         B = U7 + U8/T[i] + U9*T[i]
#         D1000 = U1*math.exp(U2*T[i] + U3*T[i]**2)
#         try:
#             D = D1000 + C*math.log((B + P[j])/(B + 1000))
#         except:
#             D = 0
#         R[i][j] = D

# plt.figure()
# plt.plot(P, R[0], 'r', P, R[1], 'b', P, R[2], 'g', P, R[3], 'r--', P, R[4], 'b--', P, R[5], 'g--', P, R[6], 'k')
# plt.legend(T)
# plt.show()

# print(f'Relative Diffusivity D = {D}') #78.38066291198618 @P=100 T=298.15