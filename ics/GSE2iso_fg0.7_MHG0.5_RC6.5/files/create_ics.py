import arepo
import numpy as np
import os
import sys
import pickle

from numba import njit

BoxSize = 400.
shift = np.array([BoxSize/2., BoxSize/2., BoxSize/2.])

Rstart = 129.
Vvir = 129.
e = 0.5
pro = 1.
angle = -165

zsolar_MW = 10.**(-0.3)
zsolar_GSE = 10.**(-1.2)

GFM_SOLAR_METALLICITY = 0.0127
GFM_INITIAL_ABUNDANCE_HYDROGEN = 0.76
GFM_INITIAL_ABUNDANCE_HELIUM = (1.-GFM_INITIAL_ABUNDANCE_HYDROGEN)
GFM_INITIAL_ABUNDANCE_CARBON = 0
GFM_INITIAL_ABUNDANCE_NITROGEN = 0
GFM_INITIAL_ABUNDANCE_OXYGEN = 0
GFM_INITIAL_ABUNDANCE_NEON = 0
GFM_INITIAL_ABUNDANCE_MAGNESIUM = 0
GFM_INITIAL_ABUNDANCE_SILICON = 0
GFM_INITIAL_ABUNDANCE_IRON = 0
GFM_INITIAL_ABUNDANCE_OTHER = 0
GFM_INITIAL_METALLICITY = 0

GFM_SOLAR_ABUNDANCE_HYDROGEN = 0.7388
GFM_SOLAR_ABUNDANCE_HELIUM  =  (1.-GFM_SOLAR_ABUNDANCE_HYDROGEN -GFM_SOLAR_METALLICITY)
GFM_SOLAR_ABUNDANCE_CARBON  =  0.0024
GFM_SOLAR_ABUNDANCE_NITROGEN =  0.0007
GFM_SOLAR_ABUNDANCE_OXYGEN  =  0.0057
GFM_SOLAR_ABUNDANCE_NEON    =  0.0012
GFM_SOLAR_ABUNDANCE_MAGNESIUM = 0.0007
GFM_SOLAR_ABUNDANCE_SILICON =  0.0007
GFM_SOLAR_ABUNDANCE_IRON   =   0.0013
GFM_SOLAR_ABUNDANCE_OTHER  =   0

def get_initial_mass_fractions(zsolar):

    mass_fractions = np.zeros(10)

    # Hydrogen
    mass_fractions[0] = GFM_INITIAL_ABUNDANCE_HYDROGEN + zsolar * (GFM_SOLAR_ABUNDANCE_HYDROGEN - GFM_INITIAL_ABUNDANCE_HYDROGEN);

    # Helium
    mass_fractions[1] = GFM_INITIAL_ABUNDANCE_HELIUM + zsolar * (GFM_SOLAR_ABUNDANCE_HELIUM - GFM_INITIAL_ABUNDANCE_HELIUM);

    # Carbon
    mass_fractions[2] = GFM_INITIAL_ABUNDANCE_CARBON + zsolar * (GFM_SOLAR_ABUNDANCE_CARBON - GFM_INITIAL_ABUNDANCE_CARBON);

    # Nitrogen
    mass_fractions[3] = GFM_INITIAL_ABUNDANCE_NITROGEN + zsolar * (GFM_SOLAR_ABUNDANCE_NITROGEN - GFM_INITIAL_ABUNDANCE_NITROGEN);

    # Oxygen
    mass_fractions[4] = GFM_INITIAL_ABUNDANCE_OXYGEN + zsolar * (GFM_SOLAR_ABUNDANCE_OXYGEN - GFM_INITIAL_ABUNDANCE_OXYGEN);

    # Neon
    mass_fractions[5] = GFM_INITIAL_ABUNDANCE_NEON + zsolar * (GFM_SOLAR_ABUNDANCE_NEON - GFM_INITIAL_ABUNDANCE_NEON);

    # Magnesium
    mass_fractions[6] = GFM_INITIAL_ABUNDANCE_MAGNESIUM + zsolar * (GFM_SOLAR_ABUNDANCE_MAGNESIUM - GFM_INITIAL_ABUNDANCE_MAGNESIUM);

    # Silicon    
    mass_fractions[7] = GFM_INITIAL_ABUNDANCE_SILICON + zsolar * (GFM_SOLAR_ABUNDANCE_SILICON - GFM_INITIAL_ABUNDANCE_SILICON);

    # Iron
    mass_fractions[8] = GFM_INITIAL_ABUNDANCE_IRON + zsolar * (GFM_SOLAR_ABUNDANCE_IRON - GFM_INITIAL_ABUNDANCE_IRON);

    # OtherMetals
    mass_fractions[9] = GFM_INITIAL_ABUNDANCE_OTHER + zsolar * (GFM_SOLAR_ABUNDANCE_OTHER - GFM_INITIAL_ABUNDANCE_OTHER);

    return mass_fractions

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

sn_MW = arepo.Snapshot('GSE_ICs.dat')
#sn_GSE = arepo.Snapshot('GSE_ICs.dat')

# create extra part type for Sgr stars
npart = np.array( [0,  0,  0,  0,  0,  0])
masses = [0., 0., 0., 0., 0., 0.]
for i in range(6):
    npart[i]  = sn_MW.NumPart_Total[i] #+ sn_GSE.NumPart_Total[i]
    masses[i] = sn_MW.MassTable[i]

ics = arepo.ICs('ics.hdf5', npart, masses=masses)

#ics.addField("GFM_Metallicity", pres=[1,0,0,0,0,0])
#ics.addField("GFM_Metals", pres=[10,0,0,0,0,0])

print(npart)
print(masses)

# Create GSE variables
#pos_GSE_part0 = sn_GSE.part0.pos
#vel_GSE_part0 = sn_GSE.part0.vel

#pos_GSE_part1 = sn_GSE.part1.pos
#vel_GSE_part1 = sn_GSE.part1.vel

#pos_GSE_part2 = sn_GSE.part2.pos
#vel_GSE_part2 = sn_GSE.part2.vel

# Offset and rotate GSE
#pos_GSE_part0, vel_GSE_part0 = separate(pos_GSE_part0, vel_GSE_part0, Rstart, Vvir, e, pro)
#pos_GSE_part1, vel_GSE_part1 = separate(pos_GSE_part1, vel_GSE_part1, Rstart, Vvir, e, pro)
#pos_GSE_part2, vel_GSE_part2 = separate(pos_GSE_part2, vel_GSE_part2, Rstart, Vvir, e, pro)

#pos_GSE_part0, vel_GSE_part0 = rotate(pos_GSE_part0, vel_GSE_part0, angle)
#pos_GSE_part1, vel_GSE_part1 = rotate(pos_GSE_part1, vel_GSE_part1, angle)
#pos_GSE_part2, vel_GSE_part2 = rotate(pos_GSE_part2, vel_GSE_part2, angle)

# Load into ics
#ics.part0.pos[:] = np.concatenate((sn_MW.part0.pos, pos_GSE_part0))
#ics.part0.mass[:] = np.concatenate((sn_MW.part0.mass, sn_GSE.part0.mass)) 
#ics.part0.vel[:] = np.concatenate((sn_MW.part0.vel, vel_GSE_part0))
#ics.part0.id[:]  = sn_gas.part0.id
#ics.part0.u[:]   = np.concatenate((sn_MW.part0.u, sn_GSE.part0.u))

ics.part0.pos[:] = sn_MW.part0.pos + shift
ics.part0.vel[:] = sn_MW.part0.vel
ics.part0.mass[:] = sn_MW.part0.mass    
ics.part0.u[:] = sn_MW.part0.u

#ics.part1.pos[:] = np.concatenate((sn_MW.part1.pos, pos_GSE_part1))
#ics.part1.vel[:] = np.concatenate((sn_MW.part1.vel, vel_GSE_part1))
#ics.part1.id[:]  = sn.part1.id

ics.part1.pos[:] = sn_MW.part1.pos + shift
ics.part1.vel[:] = sn_MW.part1.vel

#ics.part2.pos[:] = np.concatenate((sn_MW.part2.pos, pos_GSE_part2))
#ics.part2.vel[:] = np.concatenate((sn_MW.part2.vel, vel_GSE_part2))
#ics.part2.id[:]  = sn.part2.id

ics.part2.pos[:] = sn_MW.part2.pos + shift
ics.part2.vel[:] = sn_MW.part2.vel

#ics.part3.pos[:] = sn_MW.part3.pos + shift
#ics.part3.vel[:] = sn_MW.part3.vel
#ics.part3.id[:]  = sn.part3.id

# Set metallicities
#metal_fractions_MW = get_initial_mass_fractions(zsolar_MW) 
#metal_fractions_GSE = get_initial_mass_fractions(zsolar_GSE)

#metal_frac_arr_MW = np.full((sn_MW.NumPart_Total[0], len(metal_fractions_MW)), metal_fractions_MW)
#metal_frac_arr_GSE = np.full((sn_GSE.NumPart_Total[0], len(metal_fractions_GSE)), metal_fractions_GSE)

#metal_fractions = np.concatenate((metal_frac_arr_MW, metal_frac_arr_GSE))

#ics.part0.GFM_Metals[:] = metal_fractions

#metallicity_MW = np.full(sn_MW.NumPart_Total[0], zsolar_MW)
#metallicity_GSE = np.full(sn_GSE.NumPart_Total[0], zsolar_GSE)
#metallicity = np.concatenate((metallicity_MW, metallicity_GSE))

#ics.part0.GFM_Metallicity[:] = metallicity

current_id = 0
for i in range(6):
    if npart[i] > 0:
        part = getattr(ics, 'part'+str(i))
        part.id[:] = np.arange(current_id+1, current_id+1+npart[i])
        current_id = current_id + npart[i]

print('MW ', sn_MW.NumPart_Total)
#print('GSE ', sn_GSE.NumPart_Total)

#NTYPES = len(sn_MW.NumPart_Total)
#IDs = {'MW': np.zeros((NTYPES, 2), dtype='int'),
       #'GSE': np.zeros((NTYPES, 2), dtype='int')}

#for i in range(len(sn_MW.NumPart_Total)):
#    if sn_MW.NumPart_Total[i] > 0:
#        sn_MW_part = getattr(sn_MW, 'part'+str(i))
#        ics_part = getattr(ics, 'part'+str(i))
#        id_start = ics_part.id[0]
#        id_end = ics_part.id[sn_MW.NumPart_Total[i]-1]
#
#        
#
#        IDs['MW'][i][0] = id_start
#        IDs['MW'][i][1] = id_end
#    else:
#        IDs['MW'][i][0] = -1
#        IDs['MW'][i][1] = -1
#
#    if sn_GSE.NumPart_Total[i] > 0:
#        sn_MW_part = getattr(sn_MW, 'part'+str(i))
#        sn_GSE_part = getattr(sn_MW, 'part'+str(i))
#        ics_part = getattr(ics, 'part'+str(i))
#        id_start = ics_part.id[sn_MW.NumPart_Total[i]]
#        id_end = ics_part.id[sn_MW.NumPart_Total[i]+sn_GSE.NumPart_Total[i]-1]
#
#        IDs['GSE'][i][0] = id_start
#        IDs['GSE'][i][1] = id_end
#    else:
#        IDs['GSE'][i][0] = -1
#        IDs['GSE'][i][1] = -1
#
#pickle.dump(IDs, open('IDs.p', 'wb'))
#
#print('part0 MW ID: ', ics.part0.id[0], ' to ', ics.part0.id[sn_MW.NumPart_Total[0]-1])
#print('part0 GSE ID: ', ics.part0.id[sn_MW.NumPart_Total[0]], ics.part0.id[sn_MW.NumPart_Total[0] + sn_GSE.NumPart_Total[0]-1])
#print('test: ', ics.part0.id[-1])

ics.write()

