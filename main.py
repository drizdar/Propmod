import sys
import json
import promod.misc.loader as loader
data = loader.loadData(sys.argv)

print(json.dumps(data))

