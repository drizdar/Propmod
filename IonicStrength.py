import math
#Ionic Strength
pc_wt = 25
M_NaCl = 58.44277 #g/mol
Mo_s = pc_wt/100/(M_NaCl/1000)
MoNa = Mo_s #mol
MoCl = Mo_s #mol
zNa = 1 #e
zCl = -1 #e
def IonicStrength(Ions):
    I = 0
    for i in range(0, len(Ions)):
        I += 0.5*Ions[i]['molality']*math.pow(Ions[i]['charge'],2)
    return I
Ions = [
    {'molality': MoNa, 'charge': zNa},
    {'molality': MoCl, 'charge': zCl}
]
I1 = 0.5*(MoNa*zNa**2 + MoCl*zCl**2)
I2 = IonicStrength(Ions)
print(I1,I2)
print(f'Iconic Strength I = {I1}')