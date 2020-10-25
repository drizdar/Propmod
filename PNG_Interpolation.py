import json
from json import JSONEncoder
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
x_range = 3000
y_range = 90
coordinates_bounds = [[98,40],[519,41],[95,760],[518,761]]
cb = coordinates_bounds
x_start = (cb[0][0] + cb[2][0])/2
y_start = (cb[0][1] + cb[1][1])/2
height = (cb[2][1] + cb[3][1])/2-(cb[0][1] + cb[1][1])/2
width = (cb[1][0] + cb[3][0])/2-(cb[0][0] + cb[2][0])/2
print(height, width, x_start, y_start)
raw_data = [
    [98,131],
    [120,128],
    [144,124],
    [171,117],
    [193,113],
    [219,109],
    [242,103],
    [267,101],
    [291,95],
    [315,92],
    [339,87],
    [365,83],
    [390,79],
]
calc_data = np.zeros((len(raw_data),2))
for i in range(0, len(raw_data)):
    x = (raw_data[i][0] - x_start)/width*x_range
    y = (height - raw_data[i][1] + y_start)/height*y_range
    calc_data[i][0] = x
    calc_data[i][1] = y
print(calc_data)

# fig = plt.figure()
# plt.plot(calc_data[:,0], calc_data[:,1])
# plt.show()

print(calc_data.tolist())
f = open('diffusivity.json','w')
f.write(json.dumps(calc_data.tolist()))
f.close()