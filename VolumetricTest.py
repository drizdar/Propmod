from formulas import RelativeDiffusivity
import math
import json
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
report = ""
dT25 = lambda T: T - 298.15
dT0 = lambda T: T - 273.15
T = 273.15 + 25
P = 1.01325
P0 = 1.01325
U1 = 1.0837195e3
U2 = -2.4749323e-1
U3 = 1.2442861e-3
U4 = 0
U5 = -7.7222249e-2
U6 = 3.2423439e-4
U7 = -5.7917599e-7
U8 = 3.3254437e-6
U9 = 0
U10 = -2.1451068e-5
U11 = 2.2324909e-3
U12 = -6.4950599e-8
U13 = 2.4503020e-10
U14 = 0
U15 = 1.0033371e-7
U16 = -1.2784026e-6
U17 = -4.6468063e-10
U18 = 5.7054131e-13
U19 = 0
U20 = 0
U21 = 1.3581172e-10
U22 = 0
U23 = 0
U24 = -6.8152430e-6
U25 = -2.5382945e-4
U26 = 6.2480692e-8
U27 = -1.0731284e-10
U28 = 0
U = [None,U1,U2,U3,U4,U5,U6,U7,U8,U9,U10,U11,U12,U13,U14,U15,U16,U17,U18,U19,U20,U21,U22,U23,U24,U25,U26,U27,U28]
report += r'''
\begin{equation}
\begin{aligned}
    V_{(m1)} &= U1 + U2 T + U3 T^2 + U4 T^3 + (P-P0) (U5 + U6 T + U7 T^2) + \\ & (P-P0)^2 (U8 + U9 T)
    \label{eqn:Vm}
\end{aligned}
\end{equation}
\begin{align*}
    \text{Where } P &= \text{Pressure [bar]} \\ 
    P_0 &= \text{ [1.01325 bar] reference pressure} \\ 
    U1 &= \text{ Constant [$1.0837195 \times 10^{3}$]} \\ 
    U2 &= \text{ Constant [$-2.4749323 \times 10^{-1}$]} \\
    U3 &= \text{ Constant [$1.2442861 \times 10^{-3}$]} \\
    U4 &= \text{ Constant [$0$]} \\
    U5 &= \text{ Constant [$-7.7222249 \times 10^{-2}$]} \\
    U6 &= \text{ Constant [$3.2423439 \times 10^{-4}$]} \\
    U7 &= \text{ Constant [$-5.7917599 \times 10^{-7}$]} \\
    U8 &= \text{ Constant [$3.3254437 \times 10^{-6}$]} \\
    U9 &= \text{ Constant [$0$]}
\end{align*}'''
print(report)
Vm = U1 + U2*T + U3*T**2 + U4*T**3 + (P-P0)*(U5 + U6*T + U7*T**2) + (P-P0)**2*(U8 + U9*T)
print(f'Vm = {Vm}')
BV_MX = U10 + U11/(T-227) + U12*T + U13*T**2 + U14/(680-T) + (P-P0)*(U15 + U18/(T-227) + U17*T + U18**T*2 + U19/(680-T)) + (P-P0)**2*(U20 + U21/(T-227) + U22*T + U23/(680-T))
print(f'BV_MX = {BV_MX}')
TWO_CV_MX = U24 + U25/(T-227) + U26*T + U27*T**2 + U28/(680-T)
CV_MX = TWO_CV_MX/2
print(f'2CV_MX = {TWO_CV_MX}')
# vw = 1.002947
X1 = 0.99983952
X2 = 16.945176e-3
X3 = -7.987040e-6
X4 = -46.170461e-9
X5 = 105.56302e-12
X6 = -280.54253e-15
X7 = 16.879850e-3
vw = 1/((X1 + X2*dT0(T) + X3*dT0(T)**2 + X4*dT0(T)**3 + X5*dT0(T)**4 + X6*dT0(T)**5)/(1 + X7*dT0(T))) #cm^3/g
print(f'Specific volume of water {vw}')
F0 = 1.50415
F1 = 1.3421e-2
F2 = 3.0591e-5
F3 = 1.15588e-6
F4 = -5.2393e-9
F5 = 2.6561e-11
Av = F0 + F1*dT0(T) + F2*dT0(T)**2 + F3*dT0(T)**3 + F4*dT0(T)**4 + F5*dT0(T)**5
print(f'D-H Slope Volume Av = {Av}')
A0 = 16.620
A1 = 8.7385e-2
A2 = -1.9994e-3
A3 = 1.7452e-5
A4 = -0.8023e-7
V_bar_2_inf_dil = A0 + A1*(T-298.15) + A2*(T-298.15)**2 + A3*(T-298.15)**3 + A4*(T-298.15)**4
print(f'Molar volume of NaCl at infinite dilution V_bar_°NaCl = {V_bar_2_inf_dil}')
R = 83.1446261815324 #cm^3 bar / K / mol
v = 2
m = 6.011
zm = 1
zx = 1
vm = 1
vx = 1
I = m
b = 1.2
M_W = 18.01534
Y = 10
m1 = 1000/(Y*M_W)
I1 = m1
h = math.log(1+b*I**(1/2))/(2*b)
h1 = math.log(1+b*I1**(1/2))/(2*b)
M2 = 58.4428
Phi_V = V_bar_2_inf_dil + v*abs(zm*zx)*Av*h + 2*vm*vx*R*T*(m*BV_MX + m**2*(vm*zm*CV_MX))
print(f'Apparent Molar Volume Phi_V {Phi_V}')
V = (m*Phi_V +1000*vw)/(1000 + m*M2)
print(f'Specific Volume of Solution V {V}')
V_expected = 0.857301
print(1/V, V/V_expected)
V1 = m/(1000+m*M2)*(Vm/m1 + (1000/m - M_W*Y)*vw + v*abs(zm*zx)*Av*(h-h1) + 2*vm*vx*R*T*(m*BV_MX - m1*BV_MX + vm*zm*(m**2 - m1**2)*CV_MX))
report += r'''
\begin{equation}
\begin{aligned}
    V &= \frac{m}{1000+m M_2} \bigg\{ \frac{V_{(m1)}}{m1} + \left( \frac{1000}{m} - M_W Y \right) v_w + v \lvert z_m z_x \rvert A_v (h-h_1) + \\ & 2 v_m v_x R T \Big( m B^V_{MX} - m_1 B^V_{MX} + v_m z_m (m^2 - m_1^2) C^V_{MX} \Big) \bigg\}
    \label{eqn:Vm}
\end{aligned}
\end{equation}
\begin{align*}
    \text{Where } P &= \text{ Pressure [bar]} \\ 
    m &= \text{ concentration [$mol/kg_{solvent}$]} \\ 
    m_1 &=  \text{ Reference molality [$\frac{1000g}{Y M_W} mol/kg_{solvent}$]} \\ 
    Y &= 10 \text{ Number of moles water associated with each mole NaCl} \\
    M_2 &=  \text{ Molecular weight of NaCl [$58.44277 g/mol$]} \\
    v_w &= \text{ Specific volume of water at T, P} \\
    v &= 2 \text{ Number of ions NaCl dissociates to in water} \\
    v_m &= \text{ Number of cations NaCl dissociates to in water} \\
    v_x &= \text{ Number of anions NaCl dissociates to in water} \\
    z_m &= \text{ Charge of cation} \\
    z_x &= \text{ Charge of anion} \\
    B^V_{MX} &= \text{ Pressure derivative of Pitzer's $B_{MX}$ parameter} \\
    C^V_{MX} &= \text{ Pressure derivative of Pitzer's $C_{MX}$ parameter}
\end{align*}'''
print(f'Specific Volume of Solution V {V1}')
print(V1/V_expected)
print(f'Density of Solution rho {1/V1}')

#Molal @0,10,20,25 °C and 1 bar
#0.1
0.995732
0.995998
0.997620
0.998834
#0.25
0.989259
0.989781
0.991564
0.992832
#0.5
0.978889
0.979804
0.981833
0.983185
#0.75
0.968991
0.970256
0.972505
0.973932
#1
0.959525
0.961101
0.963544
0.965038
#2
0.925426
0.927905
0.930909
0.932590
#3
0.896292
0.899262
0.902565
0.904339
#4
0.870996
0.874201
0.877643
0.879457
#5
0.848646
0.851958
0.855469
0.857301
