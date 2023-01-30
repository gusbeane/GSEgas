import numpy as np
import arepo
import sys
from tqdm import tqdm
import glob
import os
import pickle
from projection import compute_projections
import h5py as h5
from numba import njit
from movie import make_movie

from joblib import Parallel, delayed

def get_pos_vel(sn, ptypes):
    pos = []
    vel = []
    for pt in ptypes:
        if sn.NumPart_Total[pt] > 0:
            part = getattr(sn, 'part'+str(pt))
            pos.append(part.pos.value)
            vel.append(part.vel.value)
    
    pos = np.concatenate(pos)
    vel = np.concatenate(vel)
    return pos, vel

def _runner(path, name, snap, ptypes=[0, 2, 3, 4], 
            rng=[[-15, 15], [-15, 15]]):
    sn = arepo.Snapshot(path + '/output/', snap, 
                        parttype=ptypes, 
                        fields=['Coordinates', 'Velocities', 'Masses'],
                        combineFiles=True)
    
    # Get pos and vel of ptypes
    pos, vel = get_pos_vel(sn, [2, 3, 4])
    pos_gas, vel_gas = get_pos_vel(sn, [0])
    
    if 'Nbody' in name:
        COM = np.array([0., 0., 0.])
    else:
        COM = np.array([200., 200., 200.])

    # Compute projection
    Hxy_s, Hxz_s, Hxy_g, Hxz_g = compute_projections(sn, COM, rng=rng)
    # Hxy_s, Hxz_s, Hxy_g, Hxz_g = None, None, None, None

    # Grab time
    time = sn.Time.value

    # Package it all together
    output = (COM, Hxy_s, Hxz_s, Hxy_g, Hxz_g, time)
    
    return output

def run(path, name, nsnap, nproc, rng):

    out = Parallel(n_jobs=nproc) (delayed(_runner)(path, name, i, rng=rng) for i in tqdm(range(nsnap)))

    COM = np.array([out[i][0] for i in range(len(out))])
    Hxy_s = np.array([out[i][1] for i in range(len(out))])
    Hxz_s = np.array([out[i][2] for i in range(len(out))])
    Hxy_g = np.array([out[i][3] for i in range(len(out))])
    Hxz_g = np.array([out[i][4] for i in range(len(out))])
    time = np.array([out[i][5] for i in range(len(out))])

    # Output non-movie stuff
    output = {'COM': COM,
              'time': time}

    pickle.dump(output, open(name+'.pickle', 'wb'))

    # Make movies
    make_movie(Hxy_s, time, Hxy_s.shape[1], 1E-3, 1E0, 'movies/'+name+'_star_xy.mp4')
    make_movie(Hxz_s, time, Hxz_s.shape[1], 1E-3, 1E0, 'movies/'+name+'_star_xz.mp4')
    make_movie(Hxy_g, time, Hxy_g.shape[1], 1E-4, 1E-1, 'movies/'+name+'_gas_xy.mp4')
    make_movie(Hxz_g, time, Hxz_g.shape[1], 1E-4, 1E-1, 'movies/'+name+'_gas_xz.mp4')

    # Output movie stuff
    # f = h5.File('proj-'+name+'.h5', 'w')
    # f.create_dataset("Hxy_s", data=Hxy_s)
    # f.create_dataset("Hxz_s", data=Hxz_s)
    # f.create_dataset("Hxy_g", data=Hxy_g)
    # f.create_dataset("Hxz_g", data=Hxz_g)
    # f.close()

if __name__ == '__main__':
    nproc = int(sys.argv[1])

    basepath = '../../runs/'

    Nbody = 'Nbody'
    fgMW05_fgGSE05 = 'fgGSE0.5_fgMW0.5'

    rng0 = [[-50, 50], [-50, 50]]
    
    pair_list = [(Nbody, 'lvl4', rng0), # 0
                 (Nbody, 'lvl3', rng0), # 1
                 (fgMW05_fgGSE05, 'lvl4', rng0), # 2
                 (fgMW05_fgGSE05, 'lvl3', rng0), # 3
                 ]

    name_list = [           p[0] + '-' + p[1] for p in pair_list]
    path_list = [basepath + p[0] + '/' + p[1] for p in pair_list]
    rng_list  = [                        p[2] for p in pair_list]

    nsnap_list = [len(glob.glob(path+'/output/snapdir*/*.0.hdf5')) for path in path_list]

  
    i = int(sys.argv[2])
    path = path_list[i]
    name = name_list[i]
    nsnap = nsnap_list[i]
    rng = rng_list[i]

    out = run(path, name, nsnap, nproc, rng)
