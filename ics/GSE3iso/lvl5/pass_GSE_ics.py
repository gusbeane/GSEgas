import arepo
import numpy as np
import os
import sys
import pickle

from numba import njit

lvl = int(sys.argv[1])
fname_in = sys.argv[2]
fname_out = sys.argv[3]

#N_GAS = {}
#N_GAS[4] = 19478 
#N_GAS[3] = 155827 
#N_GAS = N_GAS[lvl]
N_GAS = 0

BoxSize = 975 
shift = np.array([BoxSize/2., BoxSize/2., BoxSize/2.])

zsolar_GSE_disk = 10.**(-1)
zsolar_GSE_halo = 10.**(-1.2)

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

sn_GSE = arepo.Snapshot(fname_in)

npart = np.array( [0,  0,  0,  0,  0,  0])
masses = [0., 0., 0., 0., 0., 0.]
for i in range(6):
    npart[i]  = sn_GSE.NumPart_Total[i] #+ sn_GSE.NumPart_Total[i]
    masses[i] = sn_GSE.MassTable[i]

ics = arepo.ICs(fname_out, npart, masses=masses)

ics.addField("GFM_Metallicity", pres=[1,0,0,0,0,0])
ics.addField("GFM_Metals", pres=[10,0,0,0,0,0])

print(npart)
print(masses)

ics.part0.pos[:] = sn_GSE.part0.pos + shift
ics.part0.vel[:] = sn_GSE.part0.vel
ics.part0.mass[:] = sn_GSE.part0.mass    
ics.part0.u[:] = sn_GSE.part0.u

ics.part1.pos[:] = sn_GSE.part1.pos + shift
ics.part1.vel[:] = sn_GSE.part1.vel

if npart[2] > 0:
    ics.part2.pos[:] = sn_GSE.part2.pos + shift
    ics.part2.vel[:] = sn_GSE.part2.vel

if npart[3] > 0:
    ics.part3.pos[:] = sn_GSE.part3.pos + shift
    ics.part3.vel[:] = sn_GSE.part3.vel

# Set metallicities
metal_fractions_GSE_disk = get_initial_mass_fractions(zsolar_GSE_disk) 
metal_fractions_GSE_halo = get_initial_mass_fractions(zsolar_GSE_halo) 

metal_frac_arr_GSE_disk = np.full((N_GAS, 10), metal_fractions_GSE_disk)
metal_frac_arr_GSE_halo = np.full((sn_GSE.NumPart_Total[0] - N_GAS, 10), metal_fractions_GSE_halo)

metal_fractions = np.concatenate((metal_frac_arr_GSE_disk, metal_frac_arr_GSE_halo))

ics.part0.GFM_Metals[:] = metal_fractions 

metallicity_GSE_disk = np.full(N_GAS, zsolar_GSE_disk)
metallicity_GSE_halo = np.full(sn_GSE.NumPart_Total[0] - N_GAS, zsolar_GSE_halo)
metallicity = np.concatenate((metallicity_GSE_disk, metallicity_GSE_halo))

ics.part0.GFM_Metallicity[:] = metallicity

# Set ids
current_id = 0
for i in range(6):
    if npart[i] > 0:
        part = getattr(ics, 'part'+str(i))
        part.id[:] = np.arange(current_id+1, current_id+1+npart[i])
        current_id = current_id + npart[i]

#print('MW ', sn_MW.NumPart_Total)
print('GSE ', sn_GSE.NumPart_Total)

ics.write()

