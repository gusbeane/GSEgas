import numpy as np
from numba import njit
from vortrace import vortrace as vt
import logging

def get_pos_mass(sn, ptypes):
    pos = []
    mass = []
    for pt in ptypes:
        if sn.NumPart_Total[pt] > 0:
            part = getattr(sn, 'part'+str(pt))
            pos.append(part.pos.value)
            if sn.MassTable[pt] > 0.0:
                mass.append(np.full(sn.NumPart_Total[pt], sn.MassTable[pt]))
            else:
                mass.append(part.mass.value)
    
    pos = np.concatenate(pos)
    mass = np.concatenate(mass)
    return pos, mass

def _compute_gas_projection(pos, rho, BoxSize, COM, rng, nres):
    pc = vt.ProjectionCloud(pos, rho, boundbox=[0., BoxSize, 0., BoxSize, 0., BoxSize])

    extent_xy = [[COM[0] + rng[0][0], COM[0] + rng[0][1]], [COM[1] + rng[1][0], COM[1] + rng[1][1]]]
    extent_xz = [[COM[0] + rng[0][0], COM[0] + rng[0][1]], [COM[2] + rng[1][0], COM[2] + rng[1][1]]]

    bounds = [0., BoxSize]

    logging.debug('right before pc.projection')
    dat_xy = pc.projection(extent_xy, nres, bounds, COM, proj='xy')
    logging.debug('right after pc.projection')

    dat_xz = pc.projection(extent_xz, nres, bounds, COM, proj='xz')
    
    return dat_xy, dat_xz

def compute_projections(sn, COM, nres=256, rng=[[-15, 15], [-15, 15]]):
    dx = (rng[0][1]-rng[0][0]) / nres
    dy = (rng[1][1]-rng[1][0]) / nres
    surf_area = dx * dy
    
    # First do stars
    if sn.NumPart_Total[2]+sn.NumPart_Total[3]+sn.NumPart_Total[4] > 0:
        pos, mass = get_pos_mass(sn, [2, 3, 4])
        pos = pos - COM

        Hxy_s, _, _ = np.histogram2d(pos[:,0], pos[:,1], bins=(nres, nres), 
            range=rng, weights=mass/surf_area)

        Hxz_s, _, _ = np.histogram2d(pos[:,0], pos[:,2], bins=(nres, nres), 
            range=rng, weights=mass/surf_area)
    else:
        Hxy_s = Hxz_s = np.zeros((nres, nres))

    # Now do gas
    if sn.NumPart_Total[0] > 0:
        # print(dir(sn.part0))
        pos = sn.part0.pos.value
        rho = sn.part0.Density.value
        BoxSize = sn.BoxSize

        Hxy_g, Hxz_g = _compute_gas_projection(pos, rho, BoxSize, COM, rng, nres)
    else:
        Hxy_g = Hxz_g = np.zeros((nres, nres))

    return Hxy_s, Hxz_s, Hxy_g, Hxz_g
