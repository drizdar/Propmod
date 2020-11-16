import classes as cl
import formulas as f
import json
import math
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import MembraneFlow as mf
import numpy as np

T = 298.15
P = 1.01325

RO = {
    "recovery_rate": 0.4
}
RO_water = cl.flow({
    "pc_wt": 0.0001,
    "P": P,
    "T": T,
    "flow": 128e6 #L/d
})
seawater = cl.flow({
    "pc_wt": 0.035,
    "P": P,
    "T": T,
    "flow": RO_water.data["flow"] / RO["recovery_rate"] #MLD
})
seawater.CalcOsmoticProperties()
seawater.CalcMassRate()
RO_water.CalcOsmoticProperties()
RO_water.CalcMassRate()
RO_brine = cl.flow({
    "P": P,
    "T": T,
    "mass_water": seawater.data["mass_water"] - RO_water.data["mass_water"],
    "mass_NaCl": seawater.data["mass_NaCl"] - RO_water.data["mass_NaCl"],
    "mass_total": seawater.data["mass_water"] - RO_water.data["mass_water"] + \
        seawater.data["mass_NaCl"] - RO_water.data["mass_NaCl"]
})
RO_brine.CalcPcWt()
RO_brine.CalcOsmoticProperties()
RO_brine.CalcFlow()

Qb = RO_brine.GetFlow("L/d")

fs = open('membrane_data.json', 'r')
data = json.loads(fs.read())
selection_index = 0
membrane_select = data[selection_index]
membrane_select["n_leaves"] = 7 #assumed
s = 0.1
membrane = cl.membrane({
    "properties": membrane_select,
    "dimensions": {
        "Am": 40.08, #m^3
        "Ld": 1.016, #m
        "Hc": 0.0007, #m
        "Ss": 0.0045 #m
    }
},s)
membrane.CalcElementDimensions()
V = 0.05 # average velocity in membrane #m/d
gamma_D = 0.3 #ratio of draw solution to mixed solution
gamma_F = 1-gamma_D #ratio of feed solution to mixed solution
v_d = V*gamma_D/gamma_F
v_f = V*gamma_F/gamma_D
print(f'Velocity in draw channel v_d = {v_d} m/s')
print(f'Velocity in feed channel v_f = {v_f} m/s')

#k
D = 6.25e-06 #m^2/h
[d_h] = membrane.GetDimensions(["d_h"])
pc_wt = 0.057
rho = 1.06
kd, k_LHd = f.MassTransportCoefficient(D/60/60, d_h, pc_wt, rho*1000, v_d)
kf, k_LHf = f.MassTransportCoefficient(D/60/60, d_h, pc_wt, rho*1000, v_f)
print(f'Mass tranfer coefficient draw channel kd = {k_LHd}')
print(f'Mass tranfer coefficient feed channel kf = {k_LHf}')
print(f'Mass tranfer exponential factor exp(-Jw/kd) {math.exp(-6.5/k_LHd)}')
print(f'Mass tranfer exponential factor exp(-Jw/kf) {math.exp(-6.5/k_LHf)}')


#pressure loss
dPd = f.PressureLoss(membrane, pc_wt, rho, "draw", v_f)
dPf = f.PressureLoss(membrane, pc_wt, rho, "feed", v_f)
print(f'Pressure loss over membrane element draw side dPd = {dPd} bar')
print(f'Pressure loss over membrane element feed side dPf = {dPf} bar')

#calculate n_pv, Qd, mixing_ratio, and back check velocity results
[dAd, dAf, Hc, ned, nef] = membrane.GetDimensions(["dAd", "dAf", "Hc", "ned", "nef"])
[n_leaves] = membrane.GetProperties(["n_leaves"])
print(f'Number of leafs n_leaves = {n_leaves}')
print(f'Area of draw channel in membrane element dAd = {dAd} m^2')
print(f'Area of draw channel in membrane element dAf = {dAf} m^2')
print(f'Number of elements in draw flow direction per leaf ned = {ned}')
print(f'Number of elements in feed flow direction per leaf nef = {nef}')
print(f'Height of channel Hc = {Hc} m')
Qf = Qb #L/d
Qf = Qf/24/60/60/1000 #m^3/s
print(f'Total feed flow Qf = {Qf} m^3/s')
Qd_pv = v_d*dAd*nef*n_leaves #m^2/s
Qf_pv = v_f*dAf*ned #m^2/s
print(f'Draw flow per pressure vessel Qd_pv = {Qd_pv} m^3/s')
print(f'Feed flow per pressure vessel Qf_pv = {Qf_pv} m^3/s')
n_pv = math.ceil(Qf/Qf_pv)
membrane.data["properties"]["n_pv"] = n_pv
print(f'Number of pressure vessels required n_pv =  {n_pv}')
Qd = Qd_pv*n_pv
Qf = Qf_pv*n_pv
print(f'Total draw flow Qd = {Qd} m^3/s')
print(f'Total feed flow Qf = {Qf} m^3/s')
Qd_e = Qd_pv/nef/n_leaves
Qf_e = Qf_pv/ned
print(f'Draw flow per membrane element Qd_e = {Qd_e} m^3/s')
print(f'Feed flow per membrane element Qf_e = {Qf_e} m^3/s')
v_dd = Qd_e/dAd
v_ff = Qf_e/dAf
print(f'Velocity in draw channel v_d = {v_dd} m/s')
print(f'Velocity in feed channel v_f = {v_ff} m/s')
mix_ratio = Qd/Qf
print(f'Mixing ratio draw/feed flow mix_ratio = {mix_ratio}')

#L/h version
Qd_pv *= 60*60*1000
Qf_pv *= 60*60*1000
print(f'Draw flow per pressure vessel Qd_pv = {Qd_pv} L/h')
print(f'Feed flow per pressure vessel Qf_pv = {Qf_pv} L/h')
Qd = Qd_pv*n_pv
Qf = Qf_pv*n_pv
print(f'Total draw flow Qd = {Qd} L/h')
print(f'Total feed flow Qf = {Qf} L/h')
Qd_e = Qd_pv/nef/n_leaves
Qf_e = Qf_pv/ned
print(f'Draw flow per membrane element Qd_e = {Qd_e} L/h')
print(f'Feed flow per membrane element Qf_e = {Qf_e} L/h')
v_dd = Qd_e/dAd/60/60/1000
v_ff = Qf_e/dAf/60/60/1000
print(f'Velocity in draw channel v_d = {v_dd} m/s')
print(f'Velocity in feed channel v_f = {v_ff} m/s')

concentrate = cl.flow({
    "pc_wt": 0.24,
    "P": P,
    "T": T,
    "flow": Qd
})
concentrate.CalcOsmoticProperties()
concentrate.CalcMassRate()

Cd, Cf, Jw, Jw_avg, Js, Js_avg, PPd, PPf, PDens, SE = mf.calculate(membrane, concentrate, RO_brine)
print(f'Average water flux Jw_avg = {Jw_avg} L/m^2/h')
print(f'Average salt flux Js_avg = {Js_avg} g/m^2/h')
print(f'Power density of the membrane Pd = {PDens} W/m^2')
print(f'Power density of the membrane Pd = {PDens} W/m^2')
print(f'Specific energy of the solution SE = {SE} J/m^3')


[dLd, dLf, Ld, Lf, ned, nef] = membrane.GetDimensions(["dLd", "dLf", "Ld", "Lf", "ned", "nef"])
[n_leaves] = membrane.GetProperties(["n_leaves"])
Lff = Lf/n_leaves
x = np.linspace(dLd/2, Ld-dLd/2, ned)
y = np.linspace(Lff-dLf/2, dLf/2, nef)

#Jw
plt.figure()
cp = plt.contour(x[:], y[:], Jw, colors='black')
plt.clabel(cp, inline=True, fontsize=6)
plt.title(r'Water flux on membrane surface of leaf $(L\ h^{-1}\ m^{-2})$')
plt.xlabel(r'Draw position $(m)$')
plt.ylabel(r'Feed position $(m)$')
plt.savefig('./../capstone/images/memb_water_flux.png')
plt.show()

cp = plt.contour(x[:], y[:], Js, colors='black')
plt.clabel(cp, inline=True, fontsize=6)
plt.title(r'Salt flux on membrane surface of leaf $(g\ h^{-1}\ m^{-2})$')
plt.xlabel(r'Draw position $(m)$')
plt.ylabel(r'Feed position $(m)$')
plt.savefig('./../capstone/images/memb_salt_flux.png')
plt.show()

cp = plt.contour(x[:], y[:], PPd, colors='black')
plt.clabel(cp, inline=True, fontsize=6)
plt.title(r'Pressure on draw side of membrane surface of leaf $(bar)$')
plt.xlabel(r'Draw position $(m)$')
plt.ylabel(r'Feed position $(m)$')
plt.savefig('./../capstone/images/memb_pressure_draw_flux.png')
plt.show()

cp = plt.contour(x[:], y[:], PPf, colors='black')
plt.clabel(cp, inline=True, fontsize=6)
plt.title(r'Pressure on feed side of membrane surface of leaf $(bar)$')
plt.xlabel(r'Draw position $(m)$')
plt.ylabel(r'Feed position $(m)$')
plt.savefig('./../capstone/images/memb_pressure_feed_flux.png')
plt.show()

cp = plt.contour(x[:], y[:], Cd, colors='black')
plt.clabel(cp, inline=True, fontsize=6)
plt.title(r'Concentration on draw side of membrane surface of leaf $(bar)$')
plt.xlabel(r'Draw position $(m)$')
plt.ylabel(r'Feed position $(m)$')
plt.savefig('./../capstone/images/memb_conc_draw.png')
plt.show()

cp = plt.contour(x[:], y[:], Cf, colors='black')
plt.clabel(cp, inline=True, fontsize=6)
plt.title(r'Concentration on feed side of membrane surface of leaf $(bar)$')
plt.xlabel(r'Draw position $(m)$')
plt.ylabel(r'Feed position $(m)$')
plt.savefig('./../capstone/images/memb_conc_feed.png')
plt.show()