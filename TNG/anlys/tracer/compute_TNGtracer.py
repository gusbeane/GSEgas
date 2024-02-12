import numpy as np
import arepo
from tqdm import tqdm
from numba import njit
import h5py as h5
import os
import sys

from joblib import Parallel, delayed

import illustris_python as il
TNGbase = '/n/holylfs05/LABS/hernquist_lab/IllustrisTNG/Runs/L35n2160TNG/output/'

@njit
def _MC_inner_loop(ParentIDs, ParticleIDs):
    key = np.full(ParentIDs.shape, -1)
    
    ilist = np.arange(len(ParentIDs))
    jlist = np.arange(len(ParticleIDs))
    
    ilist = ilist[np.argsort(ParentIDs)]
    jlist = jlist[np.argsort(ParticleIDs)]
    
    i_ = j_ = 0
    
    while i_ < len(ilist) and j_ < len(jlist):
        if ParentIDs[ilist[i_]] == ParticleIDs[jlist[j_]]:
            key[ilist[i_]] = jlist[j_]
            i_ += 1
        elif ParentIDs[ilist[i_]] > ParticleIDs[jlist[j_]]:
            j_ += 1
        else:
            i_ += 1
    
    return key

def get_time(time, redshift=False, 
             Omega0=0.3089, 
             OmegaLambda=0.6911,
             HubbleParam=0.6774):
    HUBBLE = 3.2407789e-18
    SEC_PER_MEGAYEAR = 3.15576e13
    
    if redshift:
        a = 1./(1.+time)
    else:
        a = time
    
    fac = 2. / (3. * np.sqrt(OmegaLambda))
    ans = fac * np.arcsinh(np.sqrt(a**3 * OmegaLambda/Omega0))

    ans /= HUBBLE * HubbleParam
    ans /= SEC_PER_MEGAYEAR * 1000
    
    return ans

def _load_single_tracer_prop(chunk, fbase, fields, ParentIDs, TracerIDs):
    fname = fbase + str(chunk) + '.hdf5'
    
    fld_shape = {}
    fld_shape['Coordinates'] = 3
    fld_shape['Velocities'] = 3
    fld_shape['BirthPos'] = 3
    fld_shape['BirthVel'] = 3
    fld_shape['MagneticField'] = 3
    fld_shape['GFM_Metals'] = 10
    
    out = {}
    for fld in fields:
        if fld_shape.get(fld, 1) > 1:
            out[fld] = [np.array([]).reshape(0, fld_shape.get(fld, 1))]
        else:
            out[fld] = [np.array([])]
    
    # out['PartType'] = [np.array([]).reshape(0, 1)]
    # out['ParticleIDs'] = [np.array([]).reshape(0, 1)]
    out['PartType'] = [np.array([], dtype=np.uint8)]
    out['ParticleIDs'] = [np.array([], dtype=np.uint64)]
    out['TracerID'] = [np.array([], dtype=np.uint64)]
    
    with h5.File(fname, mode='r') as f:
        for pt in [0, 4, 5]:
            pt_str = 'PartType' + str(pt)
            
            if pt_str in f.keys():
                PtypeIDs = f[pt_str + '/ParticleIDs'][:]
                
                key0 = np.isin(PtypeIDs, ParentIDs)
                key = _MC_inner_loop(ParentIDs, PtypeIDs)
                N = len(np.where(key > -1)[0])
                N0 = len(np.where(key0)[0])
                if N0 > 0:
                    print(chunk, pt_str, N, N0)
                
                if N == 0:
                    continue
                    
                out['PartType'].append(np.full(N, pt, dtype=np.uint8))
                out['ParticleIDs'].append(PtypeIDs[key[key > -1]]) # same as ParentIDs[key > -1]
                out['TracerID'].append(TracerIDs[key > -1])
                
                for fld in fields:
                    if fld in f[pt_str].keys():
                        out[fld].append( f[pt_str + '/' + fld][:][key[key > -1]] )
                    else:
                        shape = (N, fld_shape.get(fld, 1))
                        out[fld].append( np.full(shape, np.nan) )
    
    for key in out.keys():
        out[key] = np.concatenate(out[key])
    
    return out

def load_tracer_prop(TNGbase, snapnum, fields, ParentIDs, TracerIDs, nproc=64):
    snap_str = str(snapnum).zfill(3)
    fbase = TNGbase + '/snapdir_'+ snap_str + '/snap_' + snap_str + '.'
    
    fname0 = fbase + str(0) + '.hdf5'
    with h5.File(fname0, mode='r') as f0:
        NumFiles = f0['Header'].attrs['NumFilesPerSnapshot']
    
    out = Parallel(n_jobs=nproc, backend='multiprocessing') (delayed(_load_single_tracer_prop)(i, fbase, fields, ParentIDs, TracerIDs) 
                                  for i in tqdm(range(NumFiles), position=0, leave=True))
    props = {}
    for key in out[0].keys():
        props[key] = np.concatenate([out[i][key] for i in range(len(out))])
        
    keysort = np.argsort(props['TracerID'])
    for key in props.keys():
        props[key] = props[key][keysort]
    
    return props
    
def _load_single_tracer_file(chunk, fbase, TargetIDs, match_parent=True):
    fname = fbase + str(chunk) + '.hdf5'
    
    with h5.File(fname, mode='r') as f:
        if 'PartType3' in f.keys():
            ParentID = f['PartType3/ParentID'][:]
            TracerID = f['PartType3/TracerID'][:]
            
            if match_parent:
                key = np.isin(ParentID, TargetIDs)
            else:
                key = np.isin(TracerID, TargetIDs)
            
            # key = np.isin(ParentID, TargetIDs)
            
            out = (TracerID[key], ParentID[key])
            
        else:
            out = (np.empty(0), np.empty(0))
    
    return out

def load_tracers(TNGbase, snapnum, TargetIDs, nproc=1, match_parent=True):
    snap_str = str(snapnum).zfill(3)
    fbase = TNGbase + '/snapdir_'+ snap_str + '/snap_' + snap_str + '.'
    
    fname0 = fbase + str(0) + '.hdf5'
    with h5.File(fname0, mode='r') as f0:
        NumFiles = f0['Header'].attrs['NumFilesPerSnapshot']
    
    out = Parallel(n_jobs=nproc) (delayed(_load_single_tracer_file)(i, fbase, 
                                                                    TargetIDs, match_parent=match_parent) 
                                  for i in tqdm(range(NumFiles), position=0, leave=True))
    
    TracerID = np.concatenate([out[i][0] for i in range(len(out))])
    ParentID = np.concatenate([out[i][1] for i in range(len(out))])
    
    sort = np.argsort(TracerID)
    TracerID = TracerID[sort]
    ParentID = ParentID[sort]
    
    return TracerID, ParentID

def identify_tracers(subID, nproc):
    name = 'subhalo' + str(subID)
    snap = il.snapshot.loadSubhalo(TNGbase, 99, subID, 4)
    TargetIDs = snap['ParticleIDs']
    
    TracerID, ParentID = load_tracers(TNGbase, 99, TargetIDs, nproc=nproc)
    
    fout = h5.File(name + '/tracers_' + name + '.h5', mode='w')
    fout.create_dataset('TracerID', data=TracerID)
    fout.create_dataset('ParentID', data=ParentID)
    fout.close()

    return None

def run(subID, snapnum, fields, nproc):
    name = 'subhalo' + str(subID)
    sub_dir = name
    
    # load in the tracer ids previously selected
    fin = h5.File(sub_dir + '/tracers_' + name + '.h5', mode='r')
    TargetIDs = fin['TracerID'][:]
    fin.close()
    
    # match these tracer ids to the parent ids
    TracerIDs, ParentIDs = load_tracers(TNGbase, snapnum, TargetIDs, nproc=nproc, match_parent=False)
    
    # load particle properties associated with 
    props = load_tracer_prop(TNGbase, snapnum, fields, ParentIDs, TracerIDs, nproc=64)
    
    fout = h5.File(sub_dir + '/props_' + name + '_snap' + str(snapnum).zfill(3) + '.h5', 'w')
    for key in props.keys():
        fout.create_dataset(key, data=props[key])
    
    header = arepo.Snapshot(TNGbase, snapnum, onlyHeader=True)
    fout.create_group('Header')
    fout['Header'].attrs.create('ScaleFactor', header.Time)
    fout['Header'].attrs.create('Redshift', header.Redshift)
    fout['Header'].attrs.create('Time', get_time(header.Time))
    fout['Header'].attrs.create('TimeLookBack', get_time(1) - get_time(header.Time))
    
    fout.close()
    
if __name__ == '__main__':
    nproc = int(sys.argv[1])
    subID = int(sys.argv[2])
    snapnum = int(sys.argv[3])
    
    name = 'subhalo' + str(subID)
    if not os.path.exists(name):
        os.mkdir(name)
    
    fields = ['Coordinates', 'Masses', 'GFM_Metallicity',
              'GFM_StellarFormationTime']
    
    if snapnum == -1:
        identify_tracers(subID, nproc)
    else:
        run(subID, snapnum, fields, nproc)
