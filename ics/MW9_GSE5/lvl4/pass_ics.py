import arepo
import numpy as np
import os
import sys
import pickle

from numba import njit

if len(sys.argv) != 3 and len(sys.argv) != 4:
    print('usage: ', sys.argv[0], 'fname_in fname_out [passive_scalar]')
    sys.exit(1)

fname_in = sys.argv[1]
fname_out = sys.argv[2]
if len(sys.argv) == 4:
    psc = int(sys.argv[3])


sn = arepo.Snapshot(fname_in)

BoxSize = sn.BoxSize 
shift = np.array([BoxSize/2., BoxSize/2., BoxSize/2.])

npart = np.array( [0,  0,  0,  0,  0,  0])
masses = [0., 0., 0., 0., 0., 0.]
for i in range(6):
    npart[i]  = sn.NumPart_Total[i]
    masses[i] = sn.MassTable[i]

ics = arepo.ICs(fname_out, npart, masses=masses)

ics.BoxSize = BoxSize

if npart[0] > 0:
    ics.part0.pos[:] = sn.part0.pos + shift
    ics.part0.vel[:] = sn.part0.vel
    ics.part0.mass[:] = sn.part0.mass    
    ics.part0.u[:] = sn.part0.u

    # need to clip gas pos to avoid weird box edge cases
    ics.part0.pos[:] = np.clip(ics.part0.pos, 1e-6, BoxSize*(1-1e-6))

if npart[1] > 0:
    ics.part1.pos[:] = sn.part1.pos + shift
    ics.part1.vel[:] = sn.part1.vel

if npart[2] > 0:
    ics.part2.pos[:] = sn.part2.pos + shift
    ics.part2.vel[:] = sn.part2.vel

if npart[3] > 0:
    ics.part3.pos[:] = sn.part3.pos + shift
    ics.part3.vel[:] = sn.part3.vel

# Set passive scalars
if len(sys.argv) == 4:
    ics.addField('PassiveScalars', [4, 0, 0, 0, 0, 0])
    passive = np.zeros((npart[0], 4))
    passive[:,psc] = 1.0
    ics.part0.PassiveScalars[:] = passive
    print('Passive scalar set to ', psc)
else:
    print('Passive scalar not given, passive scalar not set.')

# Set ids
current_id = 0
for i in range(6):
    if npart[i] > 0:
        part = getattr(ics, 'part'+str(i))
        part.id[:] = np.arange(current_id+1, current_id+1+npart[i])
        current_id = current_id + npart[i]

ics.write()

