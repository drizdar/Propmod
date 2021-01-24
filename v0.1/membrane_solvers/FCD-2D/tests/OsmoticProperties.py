import math
import formulas as f
T = 273.15 + 25
P = 1.01325 
pc_wt = 0.26
print(f'Percent weight {pc_wt} kg NaCl / kg Total')
m_NaCl = f.mNaCl(pc_wt)
print(f'Mass of NaCl {m_NaCl} g')
Molal_NaCl = f.Molality(m_NaCl)
# Molal_NaCl = 6
# m_NaCl = Molal_NaCl*58.44277
print(f'Molality {Molal_NaCl} mol/kg')
D = f.RelativeDiffusivity(T, P)
print(f'Relative Diffusivity D = {D}')
A_phi = f.APhi(D, T)
print(f'A phi = {A_phi}')
I = f.IonicStrength(Molal_NaCl)
print(f'Iconic Strength I = {I}')
Beta_0_NaCl = f.Beta0NaCl(T, P)
print(f'Beta 0 MX = {Beta_0_NaCl}')
Beta_1_NaCl = f.Beta1NaCl(T, P)
print(f'Beta 1 MX = {Beta_1_NaCl}')
C_phi_NaCl = f.CphiNaCl(T, P)
print(f'C phi = {C_phi_NaCl}')
print(f'C = {C_phi_NaCl/2*1e3}')
gamma_MX = f.GammaPhi(A_phi,Beta_0_NaCl,Beta_1_NaCl, C_phi_NaCl, I,Molal_NaCl)
print(f'Ionic activity coefficient gamma NaCl {gamma_MX}')
phi = f.Phi(A_phi, Beta_0_NaCl, Beta_1_NaCl, C_phi_NaCl, I, Molal_NaCl)
print(f'Osmotic coefficient phi {phi}')
A_v = f.Av(T)
print(f'A v = {A_v}')
B_V_NaCl = f.BVNaCl(T,P)
print(f'B V NaCl = {B_V_NaCl}')
C_V_NaCl = f.CVNaCl(T)
print(f'C V NaCl = {C_V_NaCl}')
V_0_NaCl = f.V0NaCl(T)
print(f'V 0 NaCl = {V_0_NaCl}')
rho_w = f.DensityWater(T)  # kg/L
print(f'Density of water = {rho_w} kg/L')
V_phi_NaCl = f.VPhiNaCL(A_v, B_V_NaCl, C_V_NaCl, I, Molal_NaCl, T, V_0_NaCl)
print(f'V phi NaCl = {V_phi_NaCl}')
rho = f.ApparentDensity(Molal_NaCl, rho_w, V_phi_NaCl)
print(f'Apparent density rho = {rho} kg/L')
V = (m_NaCl + 1000)/rho
print(f'Apparent volume = {V} cm^3')
a_w = f.WaterActivity(Molal_NaCl, phi)
print(f'Water activity coefficient a w = {a_w}')
MVW = f.MolarVolumeWater(rho_w)
print(f'Molar Volume of water {MVW} cm^3/mol')
M_NaCl = f.MolalityToMolarity(m_NaCl, rho)
print(f'Molarity of NaCl {M_NaCl} mol/L')
PI_w = f.OsmoticPressurePitzer(a_w, MVW, T)
print(f'Osmotic Pressure (Pitzer) {PI_w} bar')
PI_w2 = f.OsP(M_NaCl, 2, T)
print(f'Osmotic Pressure (van\'t Hoff) {PI_w2} bar')
print(f'Difference between van\'t Hoff and Pitzer {PI_w/PI_w2}')
PI, rho1, C = f.OsmoticProperties(P, T, pc_wt)
assert PI_w == PI
assert rho == rho1
assert C == M_NaCl
