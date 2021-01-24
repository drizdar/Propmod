import classes as cl
import formulas as f
import json

fs = open('environmental_data.json','r')
data = json.loads(fs.read())
rainfall = data.get("rainfall")
rainfall_avg = data.get("rainfall_avg")
evaporation = data.get("evaporation")
evaporation_avg = data.get("evaporation_avg")
infiltration_losses = data.get("infiltration_losses")
infiltration_losses_avg = data.get("infiltration_losses_avg")


P = 1.01325
T = 273.15 + 25

brine_pc_wt = 0.057 #%wt
brine = cl.flow({
    "P": P,
    "T": T,
    "pc_wt": brine_pc_wt,
    "flow": 192e6
})
brine.CalcOsmoticProperties()
brine.CalcMassRate()
brine_water = brine.data.get("mass_water")
brine_NaCl = brine.data.get("mass_NaCl")
rho_brine = brine.data.get("density")
m_t_brine = brine.data.get("mass_total")

targ_diff_pc_wt = 0.15 #%wt
targ_pond_pc_wt = targ_diff_pc_wt + brine_pc_wt #%wt
print(targ_pond_pc_wt) #0.207

gamma_D = 0.3
gamma_F = 1-gamma_D
print(gamma_F) #0.7

Qb = brine.GetFlow("MLD")
Qr = Qb * gamma_D/gamma_F
print(Qr) #82.28571428571428

concentrate = cl.flow({
    "P": P,
    "T": T,
    "pc_wt": targ_pond_pc_wt,
    "flow": Qr*1e6
})
concentrate.CalcOsmoticProperties()
concentrate.CalcMassRate()
concentrate_water = concentrate.data.get("mass_water")
concentrate_NaCl = concentrate.data.get("mass_NaCl")
rho_conc = concentrate.data.get("density")

discharge = cl.combineFlows(brine, concentrate)
discharge_water = discharge.data.get("mass_water")
discharge_NaCl = discharge.data.get("mass_NaCl")

rho_w = f.DensityWater(T)
rain = rainfall_avg * rho_w #kg/d/m^2
evap = evaporation_avg * rho_w #kg/d/m^2
soil_loss_water = infiltration_losses_avg * rho_conc * (1-targ_pond_pc_wt) #kg/d/m^2
d_m_w = rain - evap - soil_loss_water
# seepage_NaCl = infiltration_losses_avg * rho_conc * targ_pond_pc_wt #kg/d/m^2
# d_m_NaCl = discharge_NaCl - soil_loss_water - concentrate_NaCl

A = brine_water/d_m_w
print(A)