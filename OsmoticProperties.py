import math
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
T = 273.15 + 25
P = 0
C = U4 + U5/(U6 + T)
B = U7+ U8/T + U9*T
D1000 = U1*math.exp(U2*T + U3*T**2)
D = D1000 + C*math.log((B + P)/(B + 1000))
print(f'Relative Diffusivity D = {D}') #78.38066291198618 @P=100 T=298.15
#Debye-Hückel slope for the osmotic coefficient 
NA = 6.022045E23
PI = math.pi
pw = 1
k = 1.38066E-16 #erg K-1
# k = 1.3806049E-16 #erg K-1
e = 4.803242E-10 #esu 1/2997924580 C
A_phi   = 1/3*(2*PI*NA*pw/1000)**(1/2)*(e**2/(D*k*T))**(3/2)
print(f'A phi = {A_phi}') #0.39205622046391536
#Ionic Strength
Mo_s = 4.27 #mol 25% approx. 4.27
MoNa = Mo_s #mol
MoCl = Mo_s #mol
C_Na = 1 #e
C_Cl = -1 #e
I = 0.5*(MoNa*C_Na**2 + MoCl*C_Cl**2)
print(f'Iconic Strength I = {I}')
#Beta 0
T = 273.15 + 25
P = 0 #reference pressure 1atm
z1 = -656.81518
z2 = 24.879183
z3 = -2.1552731e-5
z4 = 5.0166855e-8
z5 = 0
z6 = -4.4640952
z7 = 0.011087099
z8 = -6.4479761e-8
z9 = -2.3234032e-10
z10 = 0
z11 = -5.2194871e-6
z12 = 2.4445210e-10
z13 = 2.8527066e-13
z14 = -1.5696231
z15 = 2.2337864e-3
z16 = -6.3933891e-7
z17 = 4.5270573e-11
z18 = 5.4151933
z19 = 0
z20 = 0
z21 = 0
z22 = 119.31966
z23 = -0.48309327
z24 = 1.4068095e-3
z25 = -4.2345814
z26 = -6.1084589
z27 = 0.40743803
z28 = -6.8152430e-6
z29 = -0.075354649
z30 = 1.2609014e-4
z31 = 6.2480692e-8
z32 = 1.8994373e-8
z33 = -1.0731264e-10
z34 = 0.32136572
z35 = -2.5382945e-4
z36 = 0
z37 = 0
Beta_0_MX = z1/T + z2 + z3*P + z4*P**2 + z5*P**3 + z6*math.log(T) + (z7+z8*P + z9*P**2 + z10*P**3)*T + (z11 + z12*P + z13*P**2)*T**2 + (z14 + z15*P + z16*P**2 + z17*P**3)/(T-227) + (z18 + z19*P + z20*P**2 + z21*P**3)/(680-T)
print(f'Beta 0 MX = {Beta_0_MX}')
#Beta 1
Beta_1_MX = z22/T + z23 + z24*T + z25/(T-227)
print(f'Beta 1 MX = {Beta_1_MX}')
#C Phi
P=0
C_phi_MX = z26/T + z27 + z28*P + z29*math.log(T) + (z30 + z31*P)*T + (z32 + z33*P)*T**2 + (z34 + z35*P)/(T-227) + (z36 + z37*P)/(680-T)
print(f'C phi = {C_phi_MX}')
print(f'C = {C_phi_MX/2*1e3}')
#Activity of water
V = 2 #dissociation constant for NaCl
Vm = 1 #cations
Vx = 1 #anions
Zm = 1 #e
Zx = -1 #e
b = 1.2 #(kg mol)^(1/2)
alpha = 2 #(kg mol)^(1/2)
ln_gamma_MX = -abs(Zm*Zx)*A_phi*(math.sqrt(I)/(1 + b*math.sqrt(1)) + 2/b*math.log(1 + b*math.sqrt(I))) + Mo_s*2*Vm*Vx/V*(2*Beta_0_MX + 2*Beta_1_MX/(alpha**2*I)*(1 - (1 + alpha*math.sqrt(I) - alpha**2*I/2)*math.exp(-alpha*math.sqrt(I)))) + 3/2*Mo_s**2*(2*(Vm*Vx)**(3/2)/V*C_phi_MX)
gamma_MX = math.exp(ln_gamma_MX)
print(f'Ionic activity coefficient gamma_MX {gamma_MX}')
phi = -abs(Zm*Zx)*A_phi*math.sqrt(I)/(1 + b*math.sqrt(1)) + Mo_s*2*Vm*Vx/V*(Beta_0_MX + Beta_1_MX*math.exp(-alpha*math.sqrt(I))) +Mo_s**2*2*(Vm*Vx)**(3/2)/V*C_phi_MX + 1
#Aparent molar volume
u1 = -2.1451068e-5
u2 = 2.2324909e-3
u3 = -6.4950599e-8
u4 = 2.4503020e-10
u5 = 0
u6 = 1.0033371e-7
u7 = -1.2784026e-6
u8 = -4.6468063e-10
u9 = 5.7054131e-13
u10 = 0
u11 = 0
u12 = 1.3581172e-10
u13 = 0
u14 = 0
u15 = -6.8152430e-6
u16 = -2.5382945e-4
u17 = 6.2480692e-8
u18 = -1.0731284e-10
u19 = 0
Av = 1.50415 + 1.3421e-2*(T-273.15) + 3.0591e-5*(T-273.15)**2 + 1.15588e-6*(T-273.15)**3 -5.2393e-9*(T-273.15)**4 + 2.6561e-11*(T-273.15)**5
print(f'Av = {Av}')
P0 = 1.01325 #bar
P = 1.01325 #bar
BV_MX = u1 + u2/(T-227) + u3*T + u4*T**2 + u5/(680-T) + (P-P0)*(u6 + u7/(T-227) + u8*T + u9**T*2 + u10/(680-T)) + (P-P0)**2*(u11 + u12/(T-227) + u13*T + u14/(680-T))
print(f'BV_MX = {BV_MX}')
TWO_CV_MX = u15 + u16/(T-227) + u17*T + u18*T**2 + u19/(680-T)
CV_MX = TWO_CV_MX/2
print(f'2CV_MX = {TWO_CV_MX}')
A0 = 16.620
A1 = 8.7385e-2
A2 = -1.9994e-3
A3 = 1.7452e-5
A4 = -0.8023e-7
V_0_MX = A0 + A1*(T-298.15) + A2*(T-298.15)**2 + A3*(T-298.15)**3 + A4*(T-298.15)**4
print(f'V_0_MX = {V_0_MX}')
M_MX = 58.44277 #g/mol
rho_w = 1000 #kg/m^3
R = 8.31446261815324e-2 #kg bar / K / mol
V_phi_MX = V_0_MX + V*abs(Zm*Zx)*Av/(2*b)*math.log(1 + b*math.sqrt(I)) + 2*R*T*Vm*Vx*Mo_s*(BV_MX + Mo_s*Vm*Zm*CV_MX)
print(f'V_phi_MX = {V_phi_MX}')
rho = (1000 + Mo_s*M_MX)/(1000/rho_w + Mo_s*V_phi_MX*1e-6)
print(f'Apparent density rho = {rho}')
m_w = rho*1 - Mo_s*M_MX
M_W = 18.01528 #g/mol
Mo_w = m_w/M_W
n_w = Mo_w*NA
m_MX = Mo_s*M_MX
pc_wt = m_MX/m_w
print(f'Percent weight solute to solvent {pc_wt}')
print(f'Moles of water n_w = {n_w}')
n_s = Mo_s*NA
T_Mo = Mo_s + Mo_w
print(Mo_s, Mo_w)
a_w = math.exp(-phi*V*Mo_s/Mo_w)
print(f'Water activity coefficient a_w = {a_w}')
PI_w = -R*T/(M_W/rho)*math.log(a_w)
PI_w2 = 2*R*T*Mo_s
print(PI_w)
print(PI_w2)
# delta_G_draw = 