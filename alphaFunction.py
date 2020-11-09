import math
import numpy as np
I_arr = np.linspace(0.1,3.5,35)
alpha = 2
for I in I_arr:
    A = 2/(alpha**2*I)*(1 - (1 + alpha*math.sqrt(I) - alpha**2*I/2)*math.exp(-alpha*math.sqrt(I)))
    print(I, A)
