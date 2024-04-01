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
#N_GAS[5] = 3120 
#N_GAS[4] = 24951 
#N_GAS[3] = #tbd
#N_GAS = N_GAS[lvl]
N_GAS = 0

BulgeRadius = 2.0
BoxSize = 1290 
shift = np.array([BoxSize/2., BoxSize/2., BoxSize/2.])

zsolar_MW_disk = 10.**(-0.3)
zsolar_MW_halo = 10.**(-1.3)

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

sn_MW = arepo.Snapshot(fname_in)

npart = np.array( [0,  0,  0,  0,  0,  0])
masses = [0., 0., 0., 0., 0., 0.]
for i in range(6):
    npart[i]  = sn_MW.NumPart_Total[i] #+ sn_GSE.NumPart_Total[i]
    masses[i] = sn_MW.MassTable[i]

ics = arepo.ICs(fname_out, npart, masses=masses)

ics.addField("GFM_Metallicity", pres=[1,0,0,0,0,0])
ics.addField("GFM_Metals", pres=[10,0,0,0,0,0])

print(npart)
print(masses)

part0_vel = np.copy(sn_MW.part0.vel)
r = np.linalg.norm(sn_MW.part0.pos, axis=1)
in_bulge = r < BulgeRadius
part0_vel[in_bulge] = 0.0

part0_u = np.copy(sn_MW.part0.u)
part0_u[in_bulge] *= 0.5

ics.part0.pos[:] = sn_MW.part0.pos + shift
ics.part0.vel[:] = part0_vel 
ics.part0.mass[:] = sn_MW.part0.mass    
ics.part0.u[:] = part0_u

ics.part1.pos[:] = sn_MW.part1.pos + shift
ics.part1.vel[:] = sn_MW.part1.vel

if npart[2] > 0:
    ics.part2.pos[:] = sn_MW.part2.pos + shift
    ics.part2.vel[:] = sn_MW.part2.vel

if npart[3] > 0:
    ics.part3.pos[:] = sn_MW.part3.pos + shift
    ics.part3.vel[:] = sn_MW.part3.vel

# Set metals
metal_fractions_MW_disk = get_initial_mass_fractions(zsolar_MW_disk) 
metal_fractions_MW_halo = get_initial_mass_fractions(zsolar_MW_halo) 

metal_frac_arr_MW_disk = np.full((N_GAS, 10), metal_fractions_MW_disk)
metal_frac_arr_MW_halo = np.full((sn_MW.NumPart_Total[0] - N_GAS, 10), metal_fractions_MW_halo)

metal_fractions = np.concatenate((metal_frac_arr_MW_disk, metal_frac_arr_MW_halo))

ics.part0.GFM_Metals[:] = metal_fractions 

# Set metallicities

metallicity_MW_disk = np.full(N_GAS, zsolar_MW_disk * GFM_SOLAR_METALLICITY)
metallicity_MW_halo = np.full(sn_MW.NumPart_Total[0] - N_GAS, zsolar_MW_halo * GFM_SOLAR_METALLICITY)
metallicity = np.concatenate((metallicity_MW_disk, metallicity_MW_halo))

ics.part0.GFM_Metallicity[:] = metallicity

# Set passive scalars
ics.addField('PassiveScalars', [4, 0, 0, 0, 0, 0])
passive_disk = np.repeat(np.array([1, 0, 0, 0], dtype=float).reshape(1, -1), N_GAS, axis=0)
passive_halo = np.repeat(np.array([0, 1, 0, 0], dtype=float).reshape(1, -1), sn_MW.NumPart_Total[0] - N_GAS, axis=0)
passive = np.concatenate((passive_disk, passive_halo))

ics.part0.PassiveScalars[:] = passive

# Set ids
current_id = 0
for i in range(6):
    if npart[i] > 0:
        part = getattr(ics, 'part'+str(i))
        part.id[:] = np.arange(current_id+1, current_id+1+npart[i])
        current_id = current_id + npart[i]

print('MW ', sn_MW.NumPart_Total)

ics.write()

