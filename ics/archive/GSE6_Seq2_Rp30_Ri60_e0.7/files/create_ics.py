import arepo
import numpy as np
import os
import sys
import pickle

from numba import njit

G = 43018.7

BoxSize = 600.
shift = np.array([BoxSize/2., BoxSize/2., BoxSize/2.])

# put Sequioa on a planar circular orbit
Ri = 60
Rp = 30
e = 0.7

pos_Seq = np.array([Ri, 0., 0.])

#sn_MW = arepo.Snapshot('MW.hdf5')
sn_GSE = arepo.Snapshot('GSE.hdf5')
sn_Seq = arepo.Snapshot('Seq.hdf5')

def get_Menc(sn, Renc):
    Menc = 0.0
    for pt in [0, 1, 2, 3, 4, 5]:
        if sn.NumPart_Total[pt] > 0:
            part = getattr(sn, 'part'+str(pt))
            r = np.linalg.norm(part.pos, axis=1)
            key = r < Renc
            
            if sn.MassTable[pt] > 0:
                Menc += len(np.where(key)[0]) * sn.MassTable[pt]
            else:
                Menc += np.sum(part.mass[key])
    return Menc

Menc = get_Menc(sn_GSE, Ri)

p = (1 + e) * Rp
phi0 = - np.arccos((p / Ri - 1)/e)
phidot = np.sqrt(G * Menc * p) / Ri / Ri
rdot = e * Ri * Ri / p * np.sin(phi0) * phidot
vphi = Ri * phidot

vel_Seq = np.array([vphi, rdot, 0])

pos_Seq_part0 = sn_Seq.part0.pos + pos_Seq
vel_Seq_part0 = sn_Seq.part0.vel + vel_Seq

pos_Seq_part1 = sn_Seq.part1.pos + pos_Seq
vel_Seq_part1 = sn_Seq.part1.vel + vel_Seq

pos_Seq_part2 = sn_Seq.part2.pos + pos_Seq
vel_Seq_part2 = sn_Seq.part2.vel + vel_Seq

# create ics
npart = np.array( [0,  0,  0,  0,  0,  0])
masses = [0., 0., 0., 0., 0., 0.]
for i in range(6):
    #npart[i]  = sn_MW.NumPart_Total[i] + sn_GSE.NumPart_Total[i]
    #masses[i] = sn_MW.MassTable[i]
    npart[i] = sn_GSE.NumPart_Total[i] + sn_Seq.NumPart_Total[i]
    masses[i] = sn_GSE.MassTable[i]

NumPart = {}
NumPart['MW'] = []
NumPart['GSE'] = []
NumPart['Seq'] = []
for i in range(6):
    NumPart['MW'].append(0)

    NumPart['GSE'].append(sn_GSE.NumPart_Total[i])

    NumPart['Seq'].append(sn_Seq.NumPart_Total[i])

#npart[0] = Npart0_MW + Npart0_GSE

ics = arepo.ICs('GSE_Seq.hdf5', npart, masses=masses)

#ics.addField("GFM_Metallicity", pres=[1,0,0,0,0,0])
#ics.addField("GFM_Metals", pres=[10,0,0,0,0,0])

print(npart)
print(masses)

# Load into ics
#ics.part0.pos[:] = np.concatenate((sn_MW.part0.pos[key_MW], pos_GSE_part0[key_GSE]))
#ics.part0.mass[:] = np.concatenate((sn_MW.part0.mass[key_MW], sn_GSE.part0.mass[key_GSE])) 
#ics.part0.vel[:] = np.concatenate((sn_MW.part0.vel[key_MW], vel_GSE_part0[key_GSE]))
#ics.part0.u[:]   = np.concatenate((sn_MW.part0.u[key_MW], sn_GSE.part0.u[key_GSE]))
#ics.part0.GFM_Metallicity[:] = np.concatenate((sn_MW.part0.GFM_Metallicity[key_MW], sn_GSE.part0.GFM_Metallicity[key_GSE]))
#ics.part0.GFM_Metals[:] = np.concatenate((sn_MW.part0.GFM_Metals[key_MW], sn_GSE.part0.GFM_Metals[key_GSE]))

#ics.part1.pos[:] = np.concatenate((sn_MW.part1.pos, pos_GSE_part1))
#ics.part1.vel[:] = np.concatenate((sn_MW.part1.vel, vel_GSE_part1))

#ics.part2.pos[:] = np.concatenate((sn_MW.part2.pos, pos_GSE_part2))
#ics.part2.vel[:] = np.concatenate((sn_MW.part2.vel, vel_GSE_part2))

#ics.part3.pos[:] = sn_MW.part3.pos
#ics.part3.vel[:] = sn_MW.part3.vel
#ics.part0.pos[:] = np.concatenate((sn_MW.part0.pos[key_MW], pos_GSE_part0[key_GSE]))

ics.part0.pos[:] = np.concatenate((sn_GSE.part0.pos, pos_Seq_part0))
ics.part0.mass[:] = np.concatenate((sn_GSE.part0.mass, sn_Seq.part0.mass)) 
ics.part0.vel[:] = np.concatenate((sn_GSE.part0.vel, vel_Seq_part0))
ics.part0.u[:]   = np.concatenate((sn_GSE.part0.u, sn_Seq.part0.u))
#ics.part0.GFM_Metallicity[:] = np.concatenate((sn_GSE.part0.GFM_Metallicity, sn_Seq.part0.GFM_Metallicity))
#ics.part0.GFM_Metals[:] = np.concatenate((sn_GSE.part0.GFM_Metals, sn_Seq.part0.GFM_Metals))

ics.part1.pos[:] = np.concatenate((sn_GSE.part1.pos, pos_Seq_part1))
ics.part1.vel[:] = np.concatenate((sn_GSE.part1.vel, vel_Seq_part1))

ics.part2.pos[:] = np.concatenate((sn_GSE.part2.pos, pos_Seq_part2))
ics.part2.vel[:] = np.concatenate((sn_GSE.part2.vel, vel_Seq_part2))

current_id = 0
for i in range(6):
    if npart[i] > 0:
        part = getattr(ics, 'part'+str(i))
        part.id[:] = np.arange(current_id+1, current_id+1+npart[i])
        current_id = current_id + npart[i]

#print('MW ', sn_MW.NumPart_Total)
print('GSE ', sn_GSE.NumPart_Total)
print('Seq ', sn_Seq.NumPart_Total)

NTYPES = len(sn_GSE.NumPart_Total)
IDs = {'MW': np.zeros((NTYPES, 2), dtype='int'),
       'GSE': np.zeros((NTYPES, 2), dtype='int'),
       'Seq': np.zeros((NTYPES, 2), dtype='int')}

for i in range(NTYPES):
    if NumPart['MW'][i] > 0:
        ics_part = getattr(ics, 'part'+str(i))
        id_start = ics_part.id[0]
        id_end = ics_part.id[NumPart['MW'][i]-1]

        IDs['MW'][i][0] = id_start
        IDs['MW'][i][1] = id_end
    else:
        IDs['MW'][i][0] = -1
        IDs['MW'][i][1] = -1

    if NumPart['GSE'][i] > 0:
        offset = NumPart['MW'][i]
        ics_part = getattr(ics, 'part'+str(i))
        id_start = ics_part.id[offset]
        id_end = ics_part.id[offset + NumPart['GSE'][i] - 1]

        IDs['GSE'][i][0] = id_start
        IDs['GSE'][i][1] = id_end
    else:
        IDs['GSE'][i][0] = -1
        IDs['GSE'][i][1] = -1
    
    if NumPart['Seq'][i] > 0:
        offset = NumPart['MW'][i] + NumPart['GSE'][i]
        ics_part = getattr(ics, 'part'+str(i))
        id_start = ics_part.id[offset]
        id_end = ics_part.id[offset + NumPart['Seq'][i] - 1]

        IDs['Seq'][i][0] = id_start
        IDs['Seq'][i][1] = id_end
    else:
        IDs['Seq'][i][0] = -1
        IDs['Seq'][i][1] = -1

pickle.dump(IDs, open('IDs.p', 'wb'))

print('part0 MW ID: ', ics.part0.id[0], ' to ', ics.part0.id[NumPart['MW'][0]-1])
print('part0 GSE ID: ', ics.part0.id[NumPart['MW'][0]], ics.part0.id[NumPart['MW'][0] + NumPart['GSE'][0] -1])
print('part0 Seq ID: ', ics.part0.id[NumPart['MW'][0] + NumPart['GSE'][0]], ics.part0.id[NumPart['MW'][0] + NumPart['GSE'][0] + NumPart['Seq'][0] - 1])
print('test: ', ics.part0.id[-1])

@njit
def my_mult(A, B):
    N = A.shape[0]
    M = B.shape[1]
    
    out = np.zeros(B.shape)
    
    for i in range(N):
        for j in range(M):
            out[i][j] = A[i] * B[i][j]
    
    return out

print('first', np.mean(ics.part0.pos, axis=0))

# subtract COM and COMV
COM = np.array([0., 0., 0.])
COMV = np.array([0., 0., 0.])
totmass = 0.
for i in range(6):
    if ics.NumPart_Total[i] > 0:
        part = getattr(ics, 'part'+str(i))
        
        if ics.MassTable[i] > 0:
            COM += np.sum(part.pos * ics.MassTable[i], axis=0)
            COMV += np.sum(part.vel * ics.MassTable[i], axis=0)
            totmass += ics.NumPart_Total[i] * ics.MassTable[i]
        else:
            COM += np.sum(my_mult(part.mass, part.pos), axis=0)
            COMV += np.sum(my_mult(part.mass, part.vel), axis=0)
            totmass += np.sum(part.mass)


COM /= totmass
COMV /= totmass

print('COM', COM)
print('COMV', COMV)

print('second', np.mean(ics.part0.pos, axis=0))

for i in range(6):
    if ics.NumPart_Total[i] > 0:
        part = getattr(ics, 'part'+str(i))
        part.pos[:] -= COM
        part.vel[:] -= COMV

print('third', np.mean(ics.part0.pos, axis=0))

ics.write()

