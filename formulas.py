from scipy.optimize import brentq
import math

# Constants
R = 0.083144598 #L-bar/mol-K
n = 2 # number of ions for Van Hoff's equation
NaCl = 58.442769 #g/mol
degCtoK = 273.15 #add to convert °C to K

# Formulas
def IDF(Qdi, Jw, ds):
    #Iterate to the next draw flow
    Qdf = Qdi + Jw*ds
    return Qdf
def IFF(Qfi, Jw, ds):
    #Iterate to the next feed flow
    Qff = Qfi - Jw*ds
    return Qff
def MF(Q, ns, npv):
    #Membrane flow
    Qm = Q/ns/npv
    return Qm
def IDC(Qdi, Cdi, Qdf, Js, ds):
    #Iterate to the next draw concentration
    Cdf = (Qdi*Cdi*NaCl - Js*ds)/Qdf/NaCl
    return Cdf
def IFC(Qfi, Cfi, Qff, Js, ds):
    #Iterate to the next feed concentration
    Cff = (Qfi*Cfi*NaCl + Js*ds)/Qff/NaCl
    return Cff
def OsP(C, n, R, T):
    #Osmotic pressure
    OsP = C*n*R*T
    return OsP
def OpP(PId, PIf):
    #Optimum pressure difference
    dP = (PId - PIf)/2
    return dP
def WF(A, B, D, k, S, PId, PIf, dP, Jwt):
    #Solve the water flux
    func = lambda Jwt: Jwt - A*((PId*math.exp(-Jwt/k) - PIf*math.exp(Jwt*1e-3*S/D))/(1 + B/Jwt*(math.exp(Jwt*1e-3*S/D)-math.exp(-Jwt/k)))-dP)
    Jw = brentq(func, 1e-6, 1000) #L/m^3/h
    return Jw
def SF(B, D, k, S, Cd, Cf, Jw):
    #Salt flux kg/m^3/h
    Js = B*(Cd*NaCl*math.exp(-Jw/k) - Cf*NaCl*math.exp(Jw*1e-3*S/D))/(1+B/Jw*(math.exp(Jw*1e-3*S/D) - math.exp(-Jw/k)))
    return Js
def V(Q, L, H):
    #Velocity in m/s
    Q = Q*1e-3/60/60 #L/h to m^3/s
    A = L * H #m^2
    V = Q/A
    return V
def BPD(Pi, Vi, Vf, Di, Df):
    #Pressure drop using Bernoulli's equation assuming density is equal
    Pf = Pi + (0.5*Vi**2*Di - 0.5*Vf**2*Df)/1e5 #Velocity component from Pa to bar
    return Pf
