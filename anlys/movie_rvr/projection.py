import numpy as np
from numba import njit
from vortrace import vortrace as vt
import logging
from scipy.stats import binned_statistic_2d

def get_r_vr(MC):
    pos = MC['PartType5/RotatedCoordinates'][:]
    vel = MC['PartType5/RotatedVelocities'][:]
    memb = MC['PartType5/Membership'][:]
    mass = MC['Header'].attrs['TracerMass']
    T = MC['PartType5/Temperature'][:]
    
    r = np.linalg.norm(pos, axis=1)
    R = np.linalg.norm(pos[:,:2], axis=1)
    cphi = pos[:,0]/R
    sphi = pos[:,1]/R
    ctheta = pos[:,2]/r
    stheta = np.sqrt(1 - ctheta**2)
    
    vr = stheta * cphi * vel[:,0] + stheta * sphi * vel[:,1] + ctheta * vel[:,2]
    
    return r[memb==1], vr[memb==1], T[memb==1], mass

def compute_projections(MC, rbins = np.linspace(20, 180, 180),
                        vrbins = np.linspace(-200, 150, 180)):

    dr = (rbins[-1]-rbins[0])/len(rbins)
    dvr = (vrbins[-1]-vrbins[0])/len(vrbins)

    r, vr, T, mass = get_r_vr(MC)
    
    m = np.full(len(r), mass/dr/dvr)
    
    Hrvr, x_edge, y_edge, _ = binned_statistic_2d(r, vr, m,
                                                  statistic='sum', bins=(rbins, vrbins))
    
    Hrvr_T, x_edge, y_edge, _ = binned_statistic_2d(r, vr, np.log10(T),
                                                  statistic='median', bins=(rbins, vrbins))

    return Hrvr, Hrvr_T
