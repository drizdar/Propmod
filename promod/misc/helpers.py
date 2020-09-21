'''
#import the tabular data and store it as a json file
import json
import TabularData as td
arr = {}
for item in td.membs:
    a = item
    arr[a.name] = {"name": a.id, "A": a.A, "B": a.B, "D": a.D, "k": a.k, "S": a.S, "material": a.mat, "reference": a.ref}
f = open('membranes.json', 'w')
f.write(json.dumps(arr))
'''
