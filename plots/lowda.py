import arepo
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import pickle
import time
from numba import njit

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

@njit
def rodrigues_formula(k, v, theta):
    N = v.shape[0]
    v_rot = np.zeros(np.shape(v))
    
    ctheta = np.cos(theta)
    stheta = np.sin(theta)
    
    for i in range(N):
        v_rot[i] = v[i] * ctheta + np.cross(k, v[i]) * stheta + k * (np.dot(k, v[i])) * (1-ctheta)
    
    return v_rot

def get_rotate_pos(pos, vel, COM, COMV, AngMom):
    pos = np.copy(pos)
    vel = np.copy(vel)
    
    pos -= COM
    vel -= COMV
    ang_mom = AngMom

    angmom_dir = ang_mom/np.linalg.norm(ang_mom)
    theta = np.arccos(np.dot(angmom_dir, np.array([0, 0, 1])))
    k = np.cross(ang_mom, np.array([0, 0, 1.]))
    k /= np.linalg.norm(k)

    pos_rot = rodrigues_formula(k, pos.astype(np.float64), theta)
    vel_rot = rodrigues_formula(k, vel.astype(np.float64), theta)
    
    return pos_rot, vel_rot

def load_gal(name, lvl, d_idx=10, idx_list=None, parttype=[0, 2, 4],
             fields = ['Coordinates', 'Velocities', 'Masses', 'PassiveScalars', 'ParticleIDs'],
             basepath='/n/holylfs05/LABS/hernquist_lab/Users/abeane/GSEgas'):
    t_start = time.time()
    
    gal = {}
    gal['sn_idx'] = {}
    gal['idx_list'] = []
    gal['i_from_idx'] = {}
    gal['Time'] = []
    
    if idx_list is None:
        idx = 0
    else:
        idx_list = np.array(idx_list)
        idx = idx_list[0]

    i = 0
    
    # convert ints into strings, but otherwise use strings
    if isinstance(lvl, int):
        lvl = 'lvl' + str(lvl)
    
    output_path = basepath + '/runs/' + name + '/' + lvl + '/output'
    
    if fields is not None:
        for fld in ['ParticleIDs', 'Coordinates', 'Masses', 'Velocities']:
            if fld not in fields:
                fields.append(fld)
    
    # load the snapshots
    while True:
        try:
            gal['sn_idx'][idx] = arepo.Snapshot(output_path, idx, parttype=parttype,
                                         fields = fields,
                                         combineFiles=True)
            gal['idx_list'].append(idx)
            gal['Time'].append(gal['sn_idx'][idx].Time.value)
            gal['i_from_idx'][idx] = i
            print('loaded idx=', idx, end='\r')
            
            if idx_list is None:
                idx += d_idx
                i += 1
            else:
                if idx == idx_list[-1]:
                    print('loaded up to idx=', idx, 'time=', gal['sn_idx'][idx].Time.value)
                    break
                
                i += 1
                idx = idx_list[i]
            
        except:
            print('loaded up to idx=', idx - d_idx, 'time=', gal['sn_idx'][idx-d_idx].Time.value)
            break
    
    gal['Time'] = np.array(gal['Time'])
    
    # compute the COM
    gal['MW_COM'] = []
    gal['MW_COMV'] = []
    gal['MW_AngMom'] = []
    gal['GSE_COM'] = []
    gal['GSE_COMV'] = []
    
    if 'iso' in name:
        IDmin = -2**63
        IDmax = 2**63 - 1
    else:
        ic_path = basepath + '/runs/' + name + '/' + lvl + '/ICs'
        IDs = pickle.load(open(ic_path + '/IDs.p', 'rb'))
        IDmin = IDs['MW'][2][0]
        IDmax = IDs['MW'][2][1]
        IDmin_GSE = IDs['GSE'][2][0]
        IDmax_GSE = IDs['GSE'][2][1]
    
    for i,idx in enumerate(gal['idx_list']):
        print('computing COM of idx=', idx, end='\r')
        sn = gal['sn_idx'][idx]
        in_MW = np.logical_and(sn.part2.id >= IDmin, sn.part2.id <= IDmax)
        
        MW_pos  = sn.part2.pos.value[in_MW]
        MW_vel  = sn.part2.vel.value[in_MW]
        MW_COM, MW_COMV = get_COM(MW_pos, MW_vel)
        MW_AngMom = get_AngMom(MW_pos, MW_vel, sn.MassTable[2], MW_COM, MW_COMV)
        
        gal['MW_COM'].append(MW_COM)
        gal['MW_COMV'].append(MW_COMV)
        gal['MW_AngMom'].append(MW_AngMom)
        
        if 'iso' not in name:
            in_GSE = np.logical_and(sn.part2.id >= IDmin_GSE, sn.part2.id <= IDmax_GSE)
            GSE_pos = sn.part2.pos.value[in_GSE]
            GSE_vel = sn.part2.vel.value[in_GSE]
            GSE_COM, GSE_COMV = get_COM(GSE_pos, GSE_vel)
            gal['GSE_COM'].append(GSE_COM)
            gal['GSE_COMV'].append(GSE_COMV)
    
        # now construct oriented coordinates
        for pt in range(6):
            if hasattr(sn, 'part'+str(pt)):
                part = getattr(sn, 'part'+str(pt))
                if hasattr(part, 'Coordinates'):
                    rotpos, rotvel = get_rotate_pos(part.pos.value, part.vel.value, MW_COM, MW_COMV, MW_AngMom)
                    part.RotatedCoordinates = rotpos
                    part.RotatedVelocities = rotvel
                    part.rotpos = rotpos
                    part.rotvel = rotvel
    
    for key in gal.keys():
        if 'COM' in key or 'AngMom' in key:
            gal[key] = np.array(gal[key])
    
    t_end = time.time()
    tdiff = t_end - t_start
    tread = time.strftime("%H:%M:%S", time.gmtime(tdiff))
    
    print("\ndone with", name, " took ", tread)

    return gal
        
