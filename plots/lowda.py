import arepo
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import pickle
import time
from numba import njit
from pathlib import Path

from joblib import Parallel, delayed
# from multiprocessing import Pool
from tqdm import tqdm

class EmptyClass:
    pass

class snap(object):
    def __init__(self, sn):
        # load header
        self.header = EmptyClass()
        for attr in dir(sn.header):
            if not attr.startswith('_'):
                setattr(self.header, attr, getattr(sn.header, attr))
        
        # load parameters
        self.parameters = EmptyClass()
        for attr in dir(sn.parameters):
            if not attr.startswith('_'):
                setattr(self.parameters, attr, getattr(sn.parameters, attr))
        
        # load parameters
        self.config = EmptyClass()
        for attr in dir(sn.config):
            if not attr.startswith('_'):
                setattr(self.config, attr, getattr(sn.config, attr))
        
        for pt in range(int(self.config.NTYPES)):
            ptype = 'part' + str(pt)
            if hasattr(sn, ptype):
                part = getattr(sn, ptype)
                setattr(self, ptype, EmptyClass())
                selfpart = getattr(self, ptype)
                for attr in dir(part):
                    if not attr.startswith('_') and not attr[0].islower():
                        setattr(selfpart, attr, np.array(getattr(part, attr)))

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

def get_idx_list(name, lvl, basepath, d_idx):
    idx_list = []
    output_path = basepath + '/runs/' + name + '/' + lvl + '/output/'
    
    idx = 0
    while True:
        idx_str = str(idx).zfill(3)
        path0 = Path(output_path + 'snapshot_' + idx_str + '.hdf5')
        path1 = Path(output_path + 'snapdir_' + idx_str + '/snapshot_' + idx_str + '.0.hdf5')
        
        if path0.exists() or path1.exists():
            idx_list.append(idx)
            idx += d_idx
        else:
            break
    
    return np.array(idx_list)

def get_n_T(sn):
    UnitLength = sn.parameters.UnitLength_in_cm
    UnitMass = sn.parameters.UnitMass_in_g
    UnitVelocity = sn.parameters.UnitVelocity_in_cm_per_s

    UnitTime = UnitLength / UnitVelocity
    UnitEnergy = UnitMass * UnitVelocity**2

    HYDROGEN_MASSFRAC = 0.76
    GAMMA = 5./3.
    PROTONMASS = 1.67262178e-24
    BOLTZMANN = 1.38065e-16

    InternalEnergy = sn.part0.InternalEnergy.value
    ElectronAbundance = sn.part0.ElectronAbundance
    Density = sn.part0.Density.value
    
    mu = 4 * PROTONMASS / (1 + 3 * HYDROGEN_MASSFRAC + 4 * HYDROGEN_MASSFRAC * ElectronAbundance)
    T = (GAMMA - 1.) * (InternalEnergy / BOLTZMANN) * (UnitEnergy / UnitMass) * mu

    n = Density / mu
    n *= UnitMass/UnitLength**3
    
    return n, T

def _load_snap(output_path, idx, parttype, fields):
    sn = arepo.Snapshot(output_path, idx, parttype=parttype, fields=fields, combineFiles=True)
    
    if 'iso' in output_path:
        IDmin = -2**63
        IDmax = 2**63 - 1
    else:
        ic_path = output_path + '/../ICs'
        # ic_path = basepath + '/runs/' + name + '/' + lvl + '/ICs'
        IDs = pickle.load(open(ic_path + '/IDs.p', 'rb'))
        IDmin = IDs['MW'][2][0]
        IDmax = IDs['MW'][2][1]
        IDmin_GSE = IDs['GSE'][2][0]
        IDmax_GSE = IDs['GSE'][2][1]
    
    in_MW = np.logical_and(sn.part2.id >= IDmin, sn.part2.id <= IDmax)
        
    MW_pos  = sn.part2.pos.value[in_MW]
    MW_vel  = sn.part2.vel.value[in_MW]
    MW_COM, MW_COMV = get_COM(MW_pos, MW_vel)
    MW_AngMom = get_AngMom(MW_pos, MW_vel, sn.MassTable[2], MW_COM, MW_COMV)
        
    if 'iso' not in output_path:
        in_GSE = np.logical_and(sn.part2.id >= IDmin_GSE, sn.part2.id <= IDmax_GSE)
        GSE_pos = sn.part2.pos.value[in_GSE]
        GSE_vel = sn.part2.vel.value[in_GSE]
        GSE_COM, GSE_COMV = get_COM(GSE_pos, GSE_vel)
    else:
        GSE_COM = np.array([np.nan, np.nan, np.nan])
        GSE_COMV = np.array([np.nan, np.nan, np.nan])
        
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
    
    n, T = get_n_T(sn)
    pres = np.zeros(sn.ntypes, dtype=int)
    pres[0] = 1
    sn.addField('Temperature', pres)
    sn.addField('NumberDensity', pres)
    sn.part0.Temperature[:] = T
    sn.part0.NumberDensity[:] = n
    
    # now construct out
    out = {}
    out['sn'] = snap(sn)
    out['Time'] = sn.Time.value
    out['MW_COM'] = MW_COM
    out['MW_COMV'] = MW_COMV
    out['MW_AngMom'] = MW_AngMom
    out['GSE_COM'] = GSE_COM
    out['GSE_COMV'] = GSE_COMV
    
    return out

def test(output_path, idx, parttype, fields):
    sn = arepo.Snapshot(output_path, idx, parttype=parttype, fields=fields, combineFiles=True)
    
    return 'asdf'

def worker(args):
    return _load_snap(*args)

def load_gal(name, lvl, d_idx=10, idx_list=None, parttype=[0, 2, 4], nproc=1,
             fields = ['Coordinates', 'Velocities', 'Masses', 'PassiveScalars', 'ParticleIDs'],
             basepath='/n/holylfs05/LABS/hernquist_lab/Users/abeane/GSEgas'):
    
    t_start = time.time()
    
    # convert ints into strings, but otherwise use strings
    if isinstance(lvl, int):
        lvl = 'lvl' + str(lvl)
    
    if idx_list is None:
        idx_list = get_idx_list(name, lvl, basepath, d_idx)
    else:
        idx_list = np.array(idx_list)

    output_path = basepath + '/runs/' + name + '/' + lvl + '/output'
    
    if fields is not None:
        for fld in ['ParticleIDs', 'Coordinates', 'Masses', 'Velocities']:
            if fld not in fields:
                fields.append(fld)
    
    outs = Parallel(n_jobs=nproc, backend='multiprocessing') (delayed(_load_snap)(output_path, idx, 
                                                      parttype, fields)
                                   for idx in tqdm(idx_list))
    
    # outs = Parallel(n_jobs=nproc) (delayed(test)(output_path, idx, 
    #                                                   parttype, fields)
    #                                for idx in tqdm(idx_list))
    
    # args = [(output_path, idx, parttype, fields) for idx in idx_list]
    
    # with Pool(nproc) as p:
        # outs = list(tqdm(p.imap(worker, args), total=len(args)))
    
    gal = {}
    gal['sn_idx'] = {}
    gal['idx_list'] = idx_list
    gal['i_from_idx'] = {}
    
    gal['Time'] = np.zeros_like(idx_list, dtype=float)
    gal['MW_COM'] = np.zeros((idx_list.shape[0], 3))
    gal['MW_COMV'] = np.zeros((idx_list.shape[0], 3))
    gal['MW_AngMom'] = np.zeros((idx_list.shape[0], 3))
    gal['GSE_COM'] = np.zeros((idx_list.shape[0], 3))
    gal['GSE_COMV'] = np.zeros((idx_list.shape[0], 3))
    
    for i,idx in enumerate(idx_list):
        out = outs[i]
        gal['sn_idx'][idx] = out['sn']
        gal['i_from_idx'][idx] = i
        
        for k in ['Time', 'MW_COM', 'MW_COMV', 'MW_AngMom', 'GSE_COM', 'GSE_COMV']:
            gal[k][i] = out[k]
    
    t_end = time.time()
    tdiff = t_end - t_start
    tread = time.strftime("%H:%M:%S", time.gmtime(tdiff))
    
    print("\ndone with", name, "got to idx=", gal['idx_list'][-1], "time=", gal['Time'][-1]) 
    print("time elapsed: ", tread)

    return gal
        
