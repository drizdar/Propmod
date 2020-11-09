import math
import formulas as f
R = 83.1446261815324  # cm^3 bar / K / mol
T = 273.15 + 25
P = 1.01325 
pc_wt = 0.26
print(f'Percent weight {pc_wt} kg NaCl / kg Total')
M_W = 1000  # mass of water g
MM_NaCl = 58.44277  # g/mol
M_NaCl = (pc_wt)/(1-pc_wt)*M_W  # mass of NaCl g
Mo_s = M_NaCl/MM_NaCl  # molality mol/kg
print(f'Molality {Mo_s} mol/kg')
MM_W = 18.01528  # g/mol
Ions = [
    {"molality": Mo_s, "charge": 1},
    {"molality": Mo_s, "charge": -1}
]
v = 2  # dissociation constant for NaCl
vm = 1  # cations
vx = 1  # anions
zm = 1  # e
zx = -1  # e
b = 1.2  # (kg / mol)^(1/2)
alpha = 2  # (kg / mol)^(1/2)
D = f.RelativeDiffusivity(T, P)
print(f'Relative Diffusivity D = {D}')
A_phi = f.APhi(D, T)
print(f'A phi = {A_phi}')
I = f.IonicStrength(Ions)
print(f'Iconic Strength I = {I}')
Beta_0_NaCl = f.Beta0NaCl(T, P)
print(f'Beta 0 MX = {Beta_0_NaCl}')
Beta_1_NaCl = f.Beta1NaCl(T, P)
print(f'Beta 1 MX = {Beta_1_NaCl}')
C_phi_NaCl = f.CphiNaCl(T, P)
print(f'C phi = {C_phi_NaCl}')
print(f'C = {C_phi_NaCl/2*1e3}')
gamma_MX = f.GammaPhi(A_phi,alpha,b, Beta_0_NaCl,Beta_1_NaCl, C_phi_NaCl, I,Mo_s,v,vm,vx,zm,zx)
print(f'Ionic activity coefficient gamma NaCl {gamma_MX}')
phi = f.phi(A_phi, alpha, b, Beta_0_NaCl, Beta_1_NaCl, C_phi_NaCl, I, Mo_s, v, vm, vx, zm, zx)
print(f'Osmotic coefficient phi {phi}')
A_v = f.Av(T)
print(f'A v = {A_v}')
B_V_NaCl = f.BVNaCl(T,P)
print(f'B V NaCl = {B_V_NaCl}')
C_V_NaCl = f.CVNaCl(T)
print(f'C V NaCl = {C_V_NaCl}')
V_0_NaCl = f.V0NaCl(T)
print(f'V 0 NaCl = {V_0_NaCl}')
rho_w = f.PW(T)  # kg/m^3 g/L
print(f'Density of water = {rho_w} kg/L')
V_phi_NaCl = f.VPhiNaCL(A_v, alpha, b, B_V_NaCl, C_V_NaCl, I, Mo_s, R, T, v, vm, vx, V_0_NaCl, zm, zx)
print(f'V phi NaCl = {V_phi_NaCl}')
rho = (1000 + Mo_s*MM_NaCl)/(1000/rho_w + Mo_s*V_phi_NaCl)
print(f'Apparent density rho = {rho} kg/L')
n_s = Mo_s
n_w = 1000/MM_W
a_w = math.exp(-phi*v*n_s/n_w)
print(f'Water activity coefficient a w = {a_w}')
VW = 1000/n_w
M_s = Mo_s/rho
print(f'Molarity of NaCl {M_s} mol/L')
PI_w = -R*T/(VW)*math.log(a_w)
print(f'Osmotic Pressure (Pitzer) {PI_w} bar')
PI_w2 = 2*R*1e-3*T*M_s
print(f'Osmotic Pressure (van\'t Hoff) {PI_w2} bar')
