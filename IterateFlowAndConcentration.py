import classes as cl
import formulas as f
import json


fs = open('membrane_data.json', 'r')
data = json.loads(fs.read())
data[0]["n_leaves"] = 7
s = 0.1
membrane = cl.membrane({
    "properties": data[0],
    "dimensions": {
        "Am": 40.08, #m^3
        "Ld": 1.016, #m
        "Hc": 0.0007, #m
        "Ss": 0.0045 #m
    }
},s)
membrane.CalcElementDimensions()
[dAd, dA] = membrane.GetDimensions(["dAd", "dAm"]) #m^2

Qi = 100 #L/h
Jw = 35 #L/h/m^2
Js = 5 #L/h/m^2
Qf = f.IDF(Qi, Jw, dA)
print(Qi, Qf)

P = 50
T = 273.15
pc_wt = 0.1
Vi = Qi/dAd/1e3/3600 #m/s

flow = cl.flow({
    "P": P,
    "T": T,
    "pc_wt": pc_wt,
    "flow": Qi

})
flow.CalcOsmoticProperties()
flow.CalcMassRate()

Ci = flow.data.get("molar_concentration")
Cf = f.IDC(Qi, Ci, Qf, Js, dA)
print(Ci, Cf)


# def IterateFlow(flow, Jw, Js, side, vel):
#     rho = flow.data.get("density")
#     rho_w = f.DensityWater(T)
#     m_w_i = flow.data.get("mass_water")
#     m_NaCl_i = flow.data.get("mass_NaCl")
#     if side == "draw":
#         m_w_f = m_w_i + Jw*dA*rho_w
#         m_NaCl_f = m_NaCl_i - Js*dA
#     else: 
#         m_w_f = m_w_i - Jw*dA*rho_w
#         m_NaCl_f = m_NaCl_i + Js*dA
#     flow.data["mass_water"] = m_w_f
#     flow.data["mass_NaCl"] = m_NaCl_f
#     flow.data["mass_total"] = m_w_f + m_NaCl_f
#     flow.CalcPcWt()
#     pc_wt = flow.data.get("pc_wt")
#     dP = f.PressureLoss(membrane, pc_wt, rho, "draw", vel)
#     flow.data["P"] -= dP
#     flow.CalcOsmoticProperties()
#     flow.CalcFlow()
#     Qdf = flow.data.get("flow")
#     Cdf = flow.data.get("molar_concentration")
#     return [Qdf, Cdf]

new_draw = f.IterateFlow(flow, membrane, Js, Jw, "draw", Vi)
Qf = new_draw.GetFlow("L/d")/24
Cf = new_draw.data["molar_concentration"]

print(Qi, Qf)
print(Ci, Cf)


Qi = 100 #L/h
Jw = 35 #L/h/m^2
Js = 5 #L/h/m^2
Qf = f.IFF(Qi, Jw, dA)
print(Qi, Qf)

P = 50
T = 273.15
pc_wt = 0.1
Vi = Qi/dAd/1e3/3600 #m/s

flow = cl.flow({
    "P": P,
    "T": T,
    "pc_wt": pc_wt,
    "flow": Qi

})
flow.CalcOsmoticProperties()
flow.CalcMassRate()

Ci = flow.data.get("molar_concentration")
Cf = f.IFC(Qi, Ci, Qf, Js, dA)
print(Ci, Cf)

new_feed = f.IterateFlow(flow, membrane, Js, Jw, "feed", Vi)
Qf = new_feed.GetFlow("L/d")/24
Cf = new_feed.data["molar_concentration"]

print(Qi, Qf)
print(Ci, Cf)

