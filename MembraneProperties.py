import formulas as f
pc_wt = 0.0578 #%wt
P = 1.01325 #bar
T = 273.15 + 25 #K
D = 6.25e-6 #m^2/h
D = D/60/60
D = 1.556e-9
PI, rho, C = f.OsmoticProperties(P, T, pc_wt) #bar, kg/L
h_c = 0.0007 #m channel height
l_f = 0.0045 #m length between spacers
Q = 12000 #L/h
A = 40.88 #m^2
L = 1.016 #m
Lf = A/L #m
vel = Q/(h_c*Lf)/3600/1000 #m/s 
d_h = f.HydraulicDiameter(h_c, l_f)
k, k_LH = f.MassTransportCoefficient(D, d_h, pc_wt, rho*1000, vel)
print(k[0], k_LH[0])
assert k == 5.72828509242054e-05
assert k_LH == 206.21826332713943