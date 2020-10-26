import math
from formulas import IonicStrength
#Ionic Strength
pc_wt = 0.25
m_total = 1 + pc_wt #kg
MM_NaCl = 58.44277 #g/mol
moNaCl = pc_wt/(MM_NaCl/1000)
moNa = moNaCl #mol
moCl = moNaCl #mol
zNa = 1 #e
zCl = -1 #e
Ions = [
    {'molality': moNa, 'charge': zNa},
    {'molality': moCl, 'charge': zCl}
]
I = IonicStrength(Ions)
assert I == 4.277689096529818, 'Ionic Strength does not match expected'
print(f'Iconic Strength I = {I}')