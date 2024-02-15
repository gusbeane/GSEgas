import arepo
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import pickle

def get_COM(pos, vel, rmin=1, rmax=20, rfac=0.9):
    COM = np.mean(pos, axis=0)
    r = np.linalg.norm(pos-COM, axis=1)
    
    rcut = rmax
    while rcut > rmin:
        COM = np.mean(pos[r<rcut], axis=0)
        r = np.linalg.norm(pos-COM, axis=1)
        rcut *= rfac
    
    COMV = np.mean(vel[r<rcut], axis=0)
    
    return COM, COMV

def get_AngMom(pos, vel, mass, COM, COMV, rcut=8):
    r = np.linalg.norm(pos-COM)
    key = r < rcut
    
    ang = np.cross(pos-COM, vel-COMV)
    ang = np.sum(ang, axis=0)
    ang *= mass
    
    return ang

def load_gal(name, lvl, d_idx=10,
             basepath='/n/holylfs05/LABS/hernquist_lab/Users/abeane/GSEgas'):
    gal = {}
    gal['sn_idx'] = {}
    gal['idx_list'] = []
    
    idx = 0
    output_path = basepath + '/runs/' + name + '/lvl' + str(lvl) + '/output'
    
    # load the snapshots
    while True:
        try:
            gal['sn_idx'][idx] = arepo.Snapshot(output_path, idx, parttype=[0, 2, 4],
                                         fields = ['Coordinates', 'Velocities', 'Masses', 
                                                   'PassiveScalars', 'ParticleIDs'],
                                         combineFiles=True)
            gal['idx_list'].append(idx)
            print('loaded idx=', idx, end='\r')
            idx += d_idx
        except:
            print('loaded up to idx=', idx - d_idx)
            break
            
    # compute the COM
    gal['MW_COM'] = {}
    
    if 'iso' in name:
        IDmin = -2**63
        IDmax = 2**63 - 1
    else:
        ic_path = basepath + '/runs/' + name + '/lvl' + str(lvl) + '/ICs'
        IDs = pickle.load(open(ic_path + '/IDs.p', 'rb'))
        IDmin = IDs['MW'][2][0]
        IDmax = IDs['MW'][2][1]
    
    for i,idx in enumerate(gal['idx_list']):
        print('computing COM of idx=', idx, end='\r')
        sn = gal['sn_idx'][idx]
        in_MW = np.logical_and(sn.part2.id >= IDmin, sn.part2.id <= IDmax)
        
        MW_pos  = sn.part2.pos.value[in_MW]
        MW_vel  = sn.part2.vel.value[in_MW]
        MW_COM, MW_COMV = get_COM(MW_pos, MW_vel)
        MW_AngMom = get_AngMom(MW_pos, MW_vel, sn.MassTable[2], MW_COM, MW_COMV)
        
        gal['MW_COM'][idx] = MW_COM
    
    print('done with', name)

    return gal
        
