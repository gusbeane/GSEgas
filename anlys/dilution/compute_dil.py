import numpy as np
import arepo
import sys
from tqdm import tqdm
import glob
import os
import pickle
import h5py as h5
from numba import njit

from joblib import Parallel, delayed

basepath = '/n/holylfs05/LABS/hernquist_lab/Users/abeane/GSEgas/'

def get_rough_metallicity(MW_mass, MW_metal_frac, GSE_mass, GSE_metal_frac):
    MW_metal = MW_mass * MW_metal_frac
    GSE_metal = GSE_mass * GSE_metal_frac
    
    return (MW_metal + GSE_metal) / (MW_mass + GSE_mass)

def get_dilution_profile(snap, name):

    fname = basepath + 'anlys/MC/'+name+'/MC_Prop_'+str(snap).zfill(3)+'.h5'
    f = h5.File(fname, mode='r')

    pos = f['PartType5/RotatedCoordinates'][:]
    ptype = f['PartType5/PartType'][:]
    mass = f['PartType5/Masses'][:]
    memb = f['PartType5/Membership'][:]
    R = np.linalg.norm(pos[:,:2], axis=1)
    z = pos[:,2]

    Rmin = 0.
    dR = 0.5

    metallicity = []
    dilution = []
    aveR = []

    while Rmin < 15:
        in_annulus = np.logical_and(R > Rmin, R < Rmin+dR)
        in_annulus = np.logical_and(in_annulus, np.abs(z) < 3)
        gas_in_annulus = np.logical_and(ptype==0, in_annulus)
    
        memb0_key = np.logical_and(gas_in_annulus, memb==0)
        memb1_key = np.logical_and(gas_in_annulus, memb==1)
        memb2_key = np.logical_and(gas_in_annulus, memb==2)
        memb3_key = np.logical_and(gas_in_annulus, memb==3)
    
        mass0 = np.sum(mass[memb0_key])
        mass1 = np.sum(mass[memb1_key])
        mass2 = np.sum(mass[memb2_key])
        mass3 = np.sum(mass[memb3_key])
    
        metal_annulus = get_rough_metallicity((mass0 + mass1), 10.**(-0.3)*0.0127,
                                              mass2+mass3, 10.**(-1.2)*0.0127)
    
        metallicity.append(metal_annulus)
        dilution.append((mass2+mass3) / (mass0 + mass1 + mass2+mass3))
        aveR.append(np.mean(R[gas_in_annulus]))

        # print(np.mean(R[gas_in_annulus]), mass0 + mass1, mass2)
        
        Rmin += dR
    
    metallicity = np.array(metallicity)
    dilution = np.array(dilution)
    aveR = np.array(aveR)
    time = f['Header'].attrs['Time']

    f.close()
    
    return time, aveR, dilution, metallicity

def _runner(name, snap, ptypes=[2]):
    
    time, aveR, dilution, metallicity = get_dilution_profile(snap, name)

    # Package it all together
    output = (time, aveR, dilution, metallicity)
    
    return output

def run(name, nsnap, nproc):

    out = Parallel(n_jobs=nproc) (delayed(_runner)(name, i) for i in tqdm(range(nsnap)))

    Time        = np.array([out[i][0] for i in range(len(out))])
    aveR        = np.array([out[i][1] for i in range(len(out))])
    dilution    = np.array([out[i][2] for i in range(len(out))])
    metallicity = np.array([out[i][3] for i in range(len(out))])

    out = {'Time'        : Time,
           'aveR'        : aveR,
           'dilution'    : dilution,
           'metallicity' : metallicity,
          }

    np.save('dil_'+name+'.npy', out)

if __name__ == '__main__':
    nproc = int(sys.argv[1])

    basepath = '../../'

    Nbody = 'Nbody'
    MW3iso_corona3 = 'MW3iso_fg0.7_MHG0.25_RC9'
    MW3_GSE2_merge2 = 'MW3_MHG0.25_GSE2_MHG0.18_Rcut10'

    pair_list = [(MW3iso_corona3, 'lvl4'), # 0
                 (MW3_GSE2_merge2, 'lvl4'), # 1
                 ]


    name_list = [           p[0] + '-' + p[1] for p in pair_list]
    path_list = [basepath + 'runs/' + p[0] + '/' + p[1] for p in pair_list]
    ic_list   = [basepath + 'ics/' + p[0] + '/' + p[1] for p in pair_list]
    
    nsnap_list = [len(glob.glob(path+'/output/snapdir*/*.0.hdf5')) for path in path_list]
  
    i = int(sys.argv[2])
    name = name_list[i]
    nsnap = nsnap_list[i]

    out = run(name, nsnap, nproc)
