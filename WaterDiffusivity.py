import math
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
#Relative Diffusivity
U1 = 3.4279E2
U2 = -5.0866E-3
U3 = 9.4690E-7
U4 = -2.0525
U5 = 3.1159E3
U6 = -1.8289E2
U7 = -8.0325E3
U8 = 4.2142E6
U9 = 2.1417
T = [298,338,373,523,573,623,673]
P = np.linspace(0,4000,200)
R = np.zeros((len(T),len(P)))
for i in range(0,len(T)):
    for j in range(0,len(P)):
        C = U4 + U5/(U6 + T[i])
        B = U7 + U8/T[i] + U9*T[i]
        D1000 = U1*math.exp(U2*T[i] + U3*T[i]**2)
        try:
            D = D1000 + C*math.log((B + P[j])/(B + 1000))
        except:
            D = 0
        R[i][j] = D

plt.figure()
plt.plot(P, R[0], 'r', P, R[1], 'b', P, R[2], 'g', P, R[3], 'r--', P, R[4], 'b--', P, R[5], 'g--', P, R[6], 'k')
plt.legend(T)
plt.show()

# print(f'Relative Diffusivity D = {D}') #78.38066291198618 @P=100 T=298.15