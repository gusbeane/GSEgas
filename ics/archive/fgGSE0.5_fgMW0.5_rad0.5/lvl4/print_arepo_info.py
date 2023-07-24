import arepo
import numpy as np
import sys

sn = arepo.Snapshot(sys.argv[1], parttype=None, fields=None, combineFiles=True)

print('MeanVolume = ', sn.BoxSize**3. / sn.NumPart_Total[0])

print(sn.NumPart_Total[0])

print((400.**3)/sn.NumPart_Total[0])

