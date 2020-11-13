import math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

import formulas as f
# pc_wt = 0.15 #%wt
pc_wt = 0.0001 #%wt
P = 1.01325
T = 273.15 +25
PI, rho, C = f.OsmoticProperties(P, T, pc_wt)
h_c = 0.0007 #m channel height
l_f = 0.0045 #m length between spacers


mu = f.InteropolateMu(pc_wt)
d_h = f.HydraulicDiameter(h_c, l_f)
vel = np.linspace(0.01, 1, 100) #m/s
Re = []
lam = []
dP = []
for i in range(0,100):
    Re.append(f.ReynoldsNumber(rho*1000, vel[i], d_h, mu))
    lam.append(f.Lambda(Re[i]))
    dP.append(f.PressureLoss(d_h, lam[i], rho, vel[i]))
print(rho, d_h, vel[0], mu, Re[0], lam[0], dP[0])



# fig, ax = plt.subplots()
# ax.plot(vel, dP, 'r')
# ax.set_xlabel('Velocity - v ($m/s$)')
# ax.set_ylabel('Pressure drop - dP ($bar/m$)')
# plt.savefig('./../capstone/images/pressure-drop.png')
# plt.show()