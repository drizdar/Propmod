import sys
import math
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import formulas as f
import MembraneElement as me
def reportInit():
    pd.set_option('display.max_rows', None)
    return np.array([["Value", "Description"]])
def pushReport(val, desc, rep):
    return np.append(rep, [[val, desc]], axis = 0)
#divide the membrane into elements
#assume one leaf per membrane
reportingOn = True
rep = reportInit()
Qd = 1800; rep = pushReport(Qd, "Draw flow into membrane (L/h)", rep) if reportingOn else rep
Qf = 500; rep = pushReport(Qf, "Feed flow into membrane (L/h)", rep) if reportingOn else rep
Am = 40.08; rep = pushReport(Am, "Area of membrane (m^2)", rep) if reportingOn else rep
Ld = 1.016; rep = pushReport(Ld, "Length of membrane along draw direction (m)", rep) if reportingOn else rep
Lf = Am/Ld; rep = pushReport(Lf, "Length of membrane along feed direction (m)", rep) if reportingOn else rep
s = 0.1; rep = pushReport(s, "Approximate length of calculation element (m)", rep) if reportingOn else rep
ned = math.ceil(Ld/s); rep = pushReport(ned, "Number of elements in draw direction", rep) if reportingOn else rep
dLd = Ld/ned; rep = pushReport(dLd, "Length of element in draw direction (m)", rep) if reportingOn else rep
nef = math.ceil(Lf/s); rep = pushReport(nef, "Number of elements in feed direction", rep) if reportingOn else rep
dLf = Lf/nef; rep = pushReport(dLf, "Length of element in feed direction (m)", rep) if reportingOn else rep
Qdi = Qd/nef; rep = pushReport(Qdi, "Draw flow along element (L/h)", rep) if reportingOn else rep
Qfi = Qf/ned; rep = pushReport(Qfi, "Feed flow along element (L/h)", rep) if reportingOn else rep
A = 0.67; rep = pushReport(A, "A: Water permeability coefficient (L/h/m^2/bar)", rep) if reportingOn else rep
B = 0.4; rep = pushReport(B, "B: Salt permeability coefficient (g/h/m^2)", rep) if reportingOn else rep
D = 6.25e-06; rep = pushReport(D, "D: Diffusion coefficient for NaCL (m^2/h)", rep) if reportingOn else rep
k = 306; rep = pushReport(k, "k: External concentration polarisation coefficient", rep) if reportingOn else rep #TODO ADD UNTIS
S = 0.0008; rep = pushReport(S, "S: Structural parameter", rep) if reportingOn else rep
n = 2
T = 25 + 273.15; rep = pushReport(T, "T: Temperature (K)", rep) if reportingOn else rep
R = f.R; rep = pushReport(R, "R: Ideal Gas Constant (mol/bar/K)", rep) if reportingOn else rep #TODO DOUBLE CHECK UNITS
Cd = np.zeros((nef, ned+1))
Pd = np.zeros((nef, ned+1))
Qd = np.zeros((nef, ned+1))
Cf = np.zeros((nef+1, ned))
Pf = np.zeros((nef+1, ned))
Qf = np.zeros((nef+1, ned))
Qd[:,0] = 18.181818181818183; rep = pushReport(Qd[0][0], "Qdi: Initial draw flow across element (L/h)", rep) if reportingOn else rep
Cd[:,0] = 3.6; rep = pushReport(Cd[0][0], "Cdi: Initial draw concentration (mol/L)", rep) if reportingOn else rep
Qf[0,:] = 47.61904761904762; rep = pushReport(Qf[0][0], "Qfi: Initial feed flow across element (L/h)", rep) if reportingOn else rep
Cf[0,:] = 0.6; rep = pushReport(Cf[0][0], "Cfi: Initial feed concentration (mol/L)", rep) if reportingOn else rep
PId = f.OsP(Cd[0][0], n, R, T); rep = pushReport(PId, "PId: Draw Osmotic Pressure (bar)", rep) if reportingOn else rep
PIf = f.OsP(Cf[0][0], n, R, T); rep = pushReport(PIf, "PIf: Feed Osmotic Pressure (bar)", rep) if reportingOn else rep
dP = (PId - PIf)/2; rep = pushReport(dP, "dP: Inital optimum pressure (bar)", rep) if reportingOn else rep
dA = dLd * dLf; rep = pushReport(dA, "dA: Area of element (m^2)", rep) if reportingOn else rep
H = 28*1e-3*0.0254; rep = pushReport(H,'H: Channel height (m)', rep) if reportingOn else rep #SW30XHR-400i assume same height for both
Pd[:,0] = 2.8 + dP; rep = pushReport(Pd[0][0],'Pdi: Initial draw pressure (bar)', rep) if reportingOn else rep
Pf[0,:] = 2.8; rep = pushReport(Pf[0][0],'Pfi: Intial feed pressure (bar)', rep) if reportingOn else rep
Jw = np.zeros((nef, ned))
Js = np.zeros((nef, ned))
Jw_pre = 10
# ret = me.calculateElement(H, dLd, dLf, A, B, D, k, S, n, R, T, Cd[0][0], Pd[0][0], Qd[0][0], Cf[0][0], Pf[0][0], Qf[0][0], 10)
# assert (ret == {'Jw': 7.7356835901890575, 'Js': 40.69055450754387, 'Cdf': 3.585574691788887, 'Pdf': 77.16868540275011, 'Qdf': 18.253175189480984, 'Cff': 0.6010355214799753, 'Pff': 2.8000007285888513, 'Qff': 47.547690611384816}), "return object"
for i in range(0, nef): #rows, along feed direction
    for j in range(0,ned): #columns, along draw direction
        ret = me.calculateElement(H, dLd, dLf, A, B, D, k, S, n, R, T, Cd[i][j], Pd[i][j], Qd[i][j], Cf[i][j], Pf[i][j], Qf[i][j], Jw_pre)
        if ret == None:
            print('exiting')
            sys.exit(1)
        else:
            Cd[i][j+1] = ret['Cdf']
            Pd[i][j+1] = ret['Pdf']
            Qd[i][j+1] = ret['Qdf']
            Cf[i+1][j] = ret['Cff']
            Pf[i+1][j] = ret['Pff']
            Qf[i+1][j] = ret['Qff']
            Jw[i][j] = ret['Jw']
            Js[i][j] = ret['Js']
            Jw_pre = ret['Jw']

print(pd.DataFrame(rep[:,0],rep[:,1],[""]))
print(np.average(Jw))
Jw_avg = np.average(Jw)
print(Cd.min())
print(Cf.max())
PDens = dP*1e5*Jw_avg/60/60/1000

print(f'Power density PD = {PDens}')
levels = MaxNLocator(nbins=50).tick_values(Jw.min(), Jw.max())


# pick the desired colormap, sensible levels, and define a normalization
# instance which takes data values and translates those into levels.
cmap = plt.get_cmap('rainbow')
norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

fig, (ax0, ax1) = plt.subplots(nrows=2)
x = np.linspace(0, Ld, ned+1)
y = np.linspace(Lf, 0, nef+1)
im = ax0.pcolormesh(x, y, Jw, cmap=cmap, norm=norm)
fig.colorbar(im, ax=ax0)
ax0.set_title('pcolormesh with levels')


# contours are *point* based plots, so convert our bound into point
# centers
cf = ax1.contourf(x[:-1] + dLd/2., y[:-1] + dLf/2., Jw, levels=levels, cmap=cmap)
fig.colorbar(cf, ax=ax1)
ax1.set_title('contourf with levels')

# adjust spacing between subplots so `ax1` title and `ax0` tick labels
# don't overlap
fig.tight_layout()

# plt.savefig('test.png')
plt.show()