import numpy as np
from numba import njit

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

def get_gas_pos_mass(sn):
    pos = sn.part0.pos.value
    mass = sn.part0.mass.value

    return pos, mass

def compute_projections(sn, COM, rng=[[-15, 15], [-15, 15]]):
    # rng = [[-15, 15], [-15, 15]]
    nres = 256

    dx = (rng[0][1]-rng[0][0]) / nres
    dy = (rng[1][1]-rng[1][0]) / nres
    surf_area = dx * dy
    
    # First do stars
    pos, mass = get_pos_mass(sn, [2, 3, 4])
    pos = pos - COM

    Hxy_s, _, _ = np.histogram2d(pos[:,0], pos[:,1], bins=(nres, nres), 
        range=rng, weights=mass/surf_area)

    Hxz_s, _, _ = np.histogram2d(pos[:,0], pos[:,2], bins=(nres, nres), 
        range=rng, weights=mass/surf_area)

    # Now do gas
    if sn.NumPart_Total[0] > 0:
        pos, mass = get_gas_pos_mass(sn)
        pos = pos - COM

        Hxy_g, _, _ = np.histogram2d(pos[:,0], pos[:,1], bins=(nres, nres), 
            range=rng, weights=mass/surf_area)

        Hxz_g, _, _ = np.histogram2d(pos[:,0], pos[:,2], bins=(nres, nres), 
            range=rng, weights=mass/surf_area)
    else:
        Hxy_g = Hxz_g = np.zeros((nres, nres))

    return Hxy_s, Hxz_s, Hxy_g, Hxz_g
