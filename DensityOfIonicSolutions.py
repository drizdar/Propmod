import math
MM_NaCl = 58.44277 #g/mol
MM_W = 18.01528 #g/mol
pc_wt = 0.342
T = 298.15
m_total = 1000 * (1 + pc_wt) #g
Mo_s = pc_wt*1000/MM_NaCl
Mo_s = 4
print(Mo_s)
MoNa = Mo_s #mol
MoCl = Mo_s #mol
n_s = Mo_s
n_w = 1000/MM_W
print(n_w)
V = 2 #dissociation constant for NaCl
VNa = 1 #cations
VCl = 1 #anions
ZNa = 1 #e
ZCl = -1 #e
b = 1.2 #(kg mol)^(1/2)
alpha = 2 #(kg mol)^(1/2)
R = 8.31446261815324 #J / K / mol
I = 0.5*(MoNa*ZNa**2 + MoCl*ZCl**2)
A0 = 16.620
A1 = 8.7385e-2
A2 = -1.9994e-3
A3 = 1.7452e-5
A4 = -0.8023e-7
dT25 = lambda T: T - 298.15
dT0 = lambda T: T - 273.15
V0NaCl = A0 + A1*dT25(T) + A2*dT25(T)**2 + A3*dT25(T)**3 + A4*dT25(T)**4
B0 = 1.2335e-5
B1 = -2.7445e-7
B2 = 2.4624e-9
B3 = -0.0108e-10
B0VNaCl = B0 + B1*dT25(T) + B2*dT25(T)**2 + B3*dT25(T)**3
C0 = 0.4354e-5
C1 = -0.9259e-6
C2 = 0.2980e-7
C3 = -0.0327e-8
B1VNaCl = C0 + C1*dT25(T) + C2*dT25(T)**2 + C3*dT25(T)**3
BVNaCl = B0VNaCl + B1VNaCl*(2/(alpha**2*I))*(1-(1+alpha*I**(1/2))*math.exp(-alpha*I**(1/2)))
E0 = -0.6578e-6
E1 = 1.5101e-8
E2 = -0.0055e-9
E3 = -0.0016e-10
CVNaCl = E0 + E1*dT25(T) + E2*dT25(T)**2 + E3*dT25(T)**3
F0 = 1.50415
F1 = 1.3421e-2
F2 = 3.0591e-5
F3 = 1.15588e-6
F4 = -5.2393e-9
F5 = 2.6561e-11
AV = F0 + F1*dT0(T) + F2*dT0(T)**2 + F3*dT0(T)**3 + F4*dT0(T)**4 + F5*dT0(T)**5
V_phi_NaCl = V0NaCl + V*abs(ZNa*ZCl)*AV/(2*b)*math.log(1 + b*math.sqrt(I)) + 2*R*T*VNa*VCl*Mo_s*(BVNaCl + Mo_s*VNa*ZNa*CVNaCl) #cm^3/mol
print(V_phi_NaCl, V0NaCl, AV, B0VNaCl, B1VNaCl, BVNaCl, CVNaCl)
P0 = (0.99983952 + 16.945176e-3*dT0(T) - 7.987040e-6*dT0(T)**2 - 46.170461e-9*dT0(T)**3 + 105.56302e-12*dT0(T)**4 - 280.54253e-15*dT0(T)**5)/(1 + 16.879850e-3*dT0(T)) #g/cm^3
print(P0,1/P0)
V_phi_NaCl = 16.62
P = (1000+Mo_s*MM_NaCl)/((1000/P0) + Mo_s*V_phi_NaCl)
print(P)
V = 1.342/P
print(V)
V_phi_NaCl = V/(n_s + n_w)
print(V_phi_NaCl)
print(0.08314*298.15*math.log(0.76)/V_phi_NaCl)