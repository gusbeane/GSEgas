import arepo
import numpy as np
import os
import sys
import pickle
from scipy.spatial import KDTree

from numba import njit

if len(sys.argv) != 4:
    print('usage: ', sys.argv[0], 'MW_ICs GSE_ICs MW_GSE_ICs')

MW_fname = sys.argv[1]
GSE_fname = sys.argv[2]
MW_GSE_fname = sys.argv[3]

sn_MW = arepo.Snapshot(MW_fname)
sn_GSE = arepo.Snapshot(GSE_fname)

BoxSize = sn_MW.BoxSize
BoxSize_GSE = sn_GSE.BoxSize

shift = np.array([BoxSize/2., BoxSize/2., BoxSize/2.])
shift_GSE = np.array([BoxSize_GSE/2., BoxSize_GSE/2., BoxSize_GSE/2.])

Rstart = 150.
Vvir = 110.
e = 0.5
pro = 1.
angle = -165

def separate(pos, vel, Rstart, Vvir, e, pro):
    vrad = np.sqrt(Vvir*Vvir - Vvir*e*Vvir*e)
    vphi = Vvir * e / np.sqrt(2.)

    pos[:,0] += Rstart

    vel[:,0] += -vrad
    vel[:,1] += pro * vphi
    vel[:,2] += vphi

    return pos, vel

def rotate(pos, vel, angle):
    theta = angle * np.pi/180.
    phi = 0.

    pos_rot = np.copy(pos)
    vel_rot = np.copy(vel)

    pos_rot[:,0] = np.cos(phi)*(np.cos(theta)*pos[:,0]+np.sin(theta)*pos[:,1])-np.sin(phi)*pos[:,1]
    pos_rot[:,1] = np.sin(phi)*(np.cos(theta)*pos[:,0]+np.sin(theta)*pos[:,2])+np.cos(phi)*pos[:,1]
    pos_rot[:,2] = -np.sin(theta)*pos[:,0]+np.cos(theta)*pos[:,2]

    vel_rot[:,0] = np.cos(phi)*(np.cos(theta)*vel[:,0]+np.sin(theta)*vel[:,1])-np.sin(phi)*vel[:,1]
    vel_rot[:,1] = np.sin(phi)*(np.cos(theta)*vel[:,0]+np.sin(theta)*vel[:,2])+np.cos(phi)*vel[:,1]
    vel_rot[:,2] = -np.sin(theta)*vel[:,0]+np.cos(theta)*vel[:,2]

    return pos_rot, vel_rot

# Create GSE variables
pos_GSE_part0 = sn_GSE.part0.pos - shift_GSE
vel_GSE_part0 = sn_GSE.part0.vel

pos_GSE_part1 = sn_GSE.part1.pos - shift_GSE
vel_GSE_part1 = sn_GSE.part1.vel

if sn_GSE.NumPart_Total[2] > 0:
    pos_GSE_part2 = sn_GSE.part2.pos - shift_GSE
    vel_GSE_part2 = sn_GSE.part2.vel

# Offset and rotate GSE
pos_GSE_part0, vel_GSE_part0 = separate(pos_GSE_part0, vel_GSE_part0, Rstart, Vvir, e, pro)
pos_GSE_part1, vel_GSE_part1 = separate(pos_GSE_part1, vel_GSE_part1, Rstart, Vvir, e, pro)
if sn_GSE.NumPart_Total[2] > 0:
    pos_GSE_part2, vel_GSE_part2 = separate(pos_GSE_part2, vel_GSE_part2, Rstart, Vvir, e, pro)

pos_GSE_part0, vel_GSE_part0 = rotate(pos_GSE_part0, vel_GSE_part0, angle)
pos_GSE_part1, vel_GSE_part1 = rotate(pos_GSE_part1, vel_GSE_part1, angle)
if sn_GSE.NumPart_Total[2] > 0:
    pos_GSE_part2, vel_GSE_part2 = rotate(pos_GSE_part2, vel_GSE_part2, angle)

pos_GSE_part0 += shift
pos_GSE_part1 += shift
if sn_GSE.NumPart_Total[2] > 0:
    pos_GSE_part2 += shift

# Get rotated pos and whether MW/GSE in Rcut_GSE
pos_cen, vel_cen = separate(np.array([[0.,0.,0.]]), np.array([[0., 0., 0.]]), Rstart, Vvir, e, pro)
pos_cen, vel_cen = rotate(pos_cen, vel_cen, angle)
pos_cen += shift

print('pos_cen =', pos_cen)

# rpart0_MW = np.linalg.norm(sn_MW.part0.pos - pos_cen, axis=1)
# rpart0_GSE = np.linalg.norm(pos_GSE_part0 - pos_cen, axis=1)

# key_MW = np.logical_not(rpart0_MW < Rcut_GSE)
# key_GSE = rpart0_GSE < Rcut_GSE

tree_MW  = KDTree(sn_MW.part0.pos)
_, key_close = tree_MW.query(pos_GSE_part0)
MW_rho  = sn_MW.part0.mass[key_close]

# key_GSE is all the ones we are keeping
key_GSE = sn_GSE.part0.mass > MW_rho
Npart0_GSE = len(np.where(key_GSE)[0])

# now, need to remove all MW cells which have their closest
# GSE cell be one that we are keeping
tree_GSE = KDTree(pos_GSE_part0)
_, key_close = tree_GSE.query(sn_MW.part0.pos)
key_MW = np.logical_not(np.isin(key_close, np.where(key_GSE)[0]))
Npart0_MW = len(np.where(key_MW)[0])

# create ics
npart = np.array( [0,  0,  0,  0,  0,  0])
masses = [0., 0., 0., 0., 0., 0.]
for i in range(6):
    npart[i]  = sn_MW.NumPart_Total[i] + sn_GSE.NumPart_Total[i]
    masses[i] = sn_MW.MassTable[i]

npart[0] = Npart0_MW + Npart0_GSE

ics = arepo.ICs(MW_GSE_fname, npart, masses=masses)

# Load into ics
ics.part0.pos[:] = np.concatenate((sn_MW.part0.pos[key_MW], pos_GSE_part0[key_GSE]))
ics.part0.mass[:] = np.concatenate((sn_MW.part0.mass[key_MW], sn_GSE.part0.mass[key_GSE])) 
ics.part0.vel[:] = np.concatenate((sn_MW.part0.vel[key_MW], vel_GSE_part0[key_GSE]))
ics.part0.u[:]   = np.concatenate((sn_MW.part0.u[key_MW], sn_GSE.part0.u[key_GSE]))

ics.part1.pos[:] = np.concatenate((sn_MW.part1.pos, pos_GSE_part1))
ics.part1.vel[:] = np.concatenate((sn_MW.part1.vel, vel_GSE_part1))

if npart[2] > 0:
    ics.part2.pos[:] = np.concatenate((sn_MW.part2.pos, pos_GSE_part2))
    ics.part2.vel[:] = np.concatenate((sn_MW.part2.vel, vel_GSE_part2))

if npart[3] > 0:
    ics.part3.pos[:] = sn_MW.part3.pos
    ics.part3.vel[:] = sn_MW.part3.vel

if hasattr(sn_MW.part0, 'PassiveScalars') and hasattr(sn_GSE.part0, 'PassiveScalars'):
    ics.addField('PassiveScalars', [4, 0, 0, 0, 0, 0])
    ics.part0.PassiveScalars[:] = np.concatenate((sn_MW.part0.PassiveScalars[key_MW], sn_GSE.part0.PassiveScalars[key_GSE]))
    print('passive scalars successfully combined')
else:
    print('missing passive scalars in one input, skipping')

current_id = 0
for i in range(6):
    if npart[i] > 0:
        part = getattr(ics, 'part'+str(i))
        part.id[:] = np.arange(current_id+1, current_id+1+npart[i])
        current_id = current_id + npart[i]

NTYPES = len(sn_MW.NumPart_Total)
IDs = {'MW': np.zeros((NTYPES, 2), dtype='int'),
       'GSE': np.zeros((NTYPES, 2), dtype='int')}

for i in range(len(sn_MW.NumPart_Total)):
    if sn_MW.NumPart_Total[i] > 0:
        if i==0:
            Npart_MW = Npart0_MW
        else:
            Npart_MW = sn_MW.NumPart_Total[i]

        sn_MW_part = getattr(sn_MW, 'part'+str(i))
        ics_part = getattr(ics, 'part'+str(i))
        id_start = ics_part.id[0]
        id_end = ics_part.id[Npart_MW-1]

        IDs['MW'][i][0] = id_start
        IDs['MW'][i][1] = id_end
    else:
        IDs['MW'][i][0] = -1
        IDs['MW'][i][1] = -1

    if sn_GSE.NumPart_Total[i] > 0:
        if i==0:
            Npart_GSE = Npart0_GSE
        else:
            Npart_GSE = sn_GSE.NumPart_Total[i]

        sn_MW_part = getattr(sn_MW, 'part'+str(i))
        sn_GSE_part = getattr(sn_MW, 'part'+str(i))
        ics_part = getattr(ics, 'part'+str(i))
        id_start = ics_part.id[Npart_MW]
        id_end = ics_part.id[Npart_MW+Npart_GSE-1]

        IDs['GSE'][i][0] = id_start
        IDs['GSE'][i][1] = id_end
    else:
        IDs['GSE'][i][0] = -1
        IDs['GSE'][i][1] = -1

pickle.dump(IDs, open('IDs.p', 'wb'))

ics.write()

