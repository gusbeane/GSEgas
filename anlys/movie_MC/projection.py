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

def get_massMC(ParticleIDs, ParentIDs, TracerMass):
    massMC = np.zeros(np.shape(ParticleIDs))

    partid_sort = np.argsort(ParticleIDs)
    parentid_sort = np.argsort(ParentIDs)

    i = j = 0
    while i < len(ParentIDs):
        parentid = ParentIDs[parentid_sort[i]]
        partid = ParticleIDs[partid_sort[j]]

        if parentid == partid:
            massMC[partid_sort[j]] += TracerMass
            i += 1
        else:
            j += 1
    
    return massMC

def compute_projections(sn, MC, COM, nres=256, rng=[[-15, 15], [-15, 15]]):
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
        # print(dir(sn.part0))
        pos = sn.part0.pos.value
        rho = sn.part0.rho.value
        mass = sn.part0.mass.value
        vol = mass/rho

        is_gas = MC['PartType5/PartType'][:] == 0
        is_GSE = MC['PartType5/Membership'][:] == 2
        is_GSE_gas = np.logical_and(is_gas, is_GSE)

        massMC = get_massMC(sn.part0.ParticleIDs, MC['PartType5/ParentID'][is_GSE_gas],
                            MC['Header'].attrs['TracerMass'])

        rhoMC = massMC / vol

        BoxSize = sn.BoxSize

        Hxy_g, Hxz_g = _compute_gas_projection(pos, rhoMC, BoxSize, COM, rng, nres)
    else:
        Hxy_g = Hxz_g = np.zeros((nres, nres))

    return Hxy_s, Hxz_s, Hxy_g, Hxz_g
