import math

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
