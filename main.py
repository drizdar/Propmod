import sys
import json
import formulas as f
from promod.misc import *
from promod.classes import *
input_data = loader.loadData(sys.argv)
file = open('promod/data/membranes.json','r')
m_data = json.loads(file.read()) #membrane data
file.close()
i_stage = input_data['stages'][0] #inital stage
n_skids = i_stage['stage']['skids'] #number of skids per stage
n_pv = i_stage['skid']['pressure_vessels'] #number of pressure vessels per skid
n_m = i_stage['pressure_vessel']['membranes'] #number of membranes per pressure vessel
t_pv = n_skids * n_pv #total number of pressure vessels
print(f'There are {t_pv} pressure vessels in this stage')
initial_draw = input_data['draw']
initial_feed = input_data['feed']
units = input_data['units']
flow_units = units['flow']
pressure_vessel_feed_flow = initial_feed['flow'] / t_pv
print(f'The flow through each pressure vessel is {pressure_vessel_feed_flow:.3} {flow_units}')
m_id = i_stage['membrane']['id']
m_prop = m_data[m_id]
m_dims = i_stage['membrane']['dimensions']
print(f'Analysis using {m_id}')
A = m_prop.get('A')
B = m_prop.get('B')
D = m_prop.get('D')
S = m_prop.get('S')
k = m_prop.get('k')
Am = m_dims.get('active_area')
Lm = m_dims.get('length')
W = Am/Lm
TAm = Am*n_m
# Initial conditions
n_eppm = 10 # number of evaluation points per membrane
dA = Am/n_eppm
t_ep = n_eppm * n_m
Td = initial_draw['temperature'] + f.degCtoK
Qd = [initial_draw['flow']/ ]
Cd = [initial_draw['concentration']]
PId = [f.OsP(Cd[0], f.n, f.R, Td)]
Tf = initial_feed['temperature'] + f.degCtoK
Qf = [initial_feed['flow']]
Cf = [initial_feed['concentration']]
PIf = [f.OsP(Cf[0], f.n, f.R, Tf)]
dP = f.OpP(PId[0], PIf[0])
Jw = [f.WF(A, B, D, k, S, PId[0], PIf[0],dP, 10)]
Js = [f.SF(B, D, k, S, Cd[0], Cf[0], Jw[0])]
for i in range(0, t_ep):
    Qd.append(f.IDF(Qd[i], Jw[i], dA))
    Qf.append(f.IFF(Qf[i], Jw[i], dA))
    Cd.append(f.IDC(Qd[i], Cd[i], Qd[i+1], Js[i], dA))
    Cf.append(f.IFC(Qf[i], Cf[i], Qf[i+1], Js[i], dA))
    PId.append(f.OsP(Cd[i+1], f.n, f.R, Td))
    PIf.append(f.OsP(Cf[i+1], f.n, f.R, Tf))
    Jw.append(f.WF(A, B, D, k, S, PId[i+1], PIf[i+1],dP, Jw[i]))
    Js.append(f.SF(B, D, k, S, Cd[i+1], Cf[i+1], Jw[i+1]))
    if (Jw[i+1]) < 0.01:
        print('Flux dropped too low')
        break
print(Jw)
print(Cd)