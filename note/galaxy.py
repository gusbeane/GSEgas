import numpy as np
import arepo
from numba import njit

@njit
def rodrigues_formula(k, v, theta):
    N = v.shape[0]
    v_rot = np.zeros(np.shape(v))
    
    ctheta = np.cos(theta)
    stheta = np.sin(theta)
    
    for i in range(N):
        v_rot[i] = v[i] * ctheta + np.cross(k, v[i]) * stheta + k * (np.dot(k, v[i])) * (1-ctheta)
    
    return v_rot

class CenterOfMass(object):
    pass

class Galaxy(object):
    def __init__(self, output_dir, snapnum, parttype=None, fields=None, orient=False, 
                 subID=0, rhalf_fac=1):
        self.sn = arepo.Snapshot(output_dir, snapnum, 
                                 parttype=parttype, fields=None, combineFiles=True)
        self.sub = arepo.Subfind(output_dir, snapnum, combineFiles=True)
        
        self.compute_n_T()
        
        if orient:
            self.do_orient(rhalf_fac=rhalf_fac)
        
    def compute_n_T(self):
        UnitLength = self.sn.parameters.UnitLength_in_cm
        UnitMass = self.sn.parameters.UnitMass_in_g
        UnitVelocity = self.sn.parameters.UnitVelocity_in_cm_per_s

        UnitTime = UnitLength / UnitVelocity
        UnitEnergy = UnitMass * UnitVelocity**2

        HYDROGEN_MASSFRAC = 0.76
        GAMMA = 5./3.
        PROTONMASS = 1.67262178e-24
        BOLTZMANN = 1.38065e-16

        InternalEnergy = self.sn.part0.InternalEnergy.value
        if hasattr(self.sn.part0, 'ElectronAbundance'):
            ElectronAbundance = self.sn.part0.ElectronAbundance
        else:
            ElectronAbundance = np.full_like(sn.part0.Density.value, (2+5*HYDROGEN_MASSFRAC)/(4*HYDROGEN_MASSFRAC))
        
        Density = self.sn.part0.Density.value
    
        mu = 4 * PROTONMASS / (1 + 3 * HYDROGEN_MASSFRAC + 4 * HYDROGEN_MASSFRAC * ElectronAbundance)
        T = (GAMMA - 1.) * (InternalEnergy / BOLTZMANN) * (UnitEnergy / UnitMass) * mu

        n = Density / mu
        n *= UnitMass/UnitLength**3
        
        self.sn.part0.HydrogenNumberDensity = np.copy(n)
        self.sn.part0.nH = self.sn.part0.HydrogenNumberDensity.view()
        
        self.sn.addField('Temperature', [1, 0, 0, 0, 0, 0])
        self.sn.part0.Temperature[:] = np.copy(T)
        self.sn.part0.T = self.sn.part0.Temperature.view()
        
    def do_orient(self, subID=0, rhalf_fac=1):
        if len(self.sub.SubhaloPos) <= subID:
            print('subID=', subID, 'doesnt exist')
            return
        
        # don't orient, just use box frame
        if self.sn.NumPart_Total[4] < 64:
            COM = self.sub.SubhaloPos[subID]
            COMV = self.sub.SubhaloVel[subID]
            
            for pt,Npart in enumerate(self.sn.NumPart_Total):
                if Npart == 0:
                    continue
            
                part = getattr(self.sn, 'part'+str(pt))
                pos = part.Coordinates.value - COM
                vel = part.Velocities.value - COMV
    
                pos_rot = np.copy(pos)
                vel_rot = np.copy(vel)
        
                part.RotatedCoordinates = np.copy(pos_rot)
                part.RotatedVelocities = np.copy(vel_rot)
        
                part.rotpos = part.RotatedCoordinates.view()
                part.rotvel = part.RotatedVelocities.view()
            
            return
        
        # get COM as subhalo pos
        COM = self.sub.SubhaloPos[subID]
        pos = self.sn.part4.pos.value - COM
        r = np.linalg.norm(pos, axis=1)
        
        # get stars within rhalf_fac * rhalf of COM
        rhalf = self.sub.SubhaloHalfmassRadType[subID,4]
        in_rhalf = r < rhalf_fac * rhalf
        is_star = self.sn.part4.GFM_StellarFormationTime > 0
        is_star_in_rhalf = np.logical_and(is_star, in_rhalf)
        
        # get COMV as mass weighted vel of stars within rhalf_fac * rhalf of COM
        vel_in_rhalf = self.sn.part4.Velocities[is_star_in_rhalf]
        mass_in_rhalf = self.sn.part4.Masses[is_star_in_rhalf]
        COMV = np.average(vel_in_rhalf, axis=0, weights=mass_in_rhalf)
    
        # compute ang mom
        vel = self.sn.part4.vel.value - COMV
        AngMom = np.cross(pos[is_star_in_rhalf], vel[is_star_in_rhalf])
        AngMom *= mass_in_rhalf.reshape(-1, 1)
        AngMom = np.sum(AngMom, axis=0)
        
        self.CenterOfMass = CenterOfMass()
        self.CenterOfMass.Coordinate = COM
        self.CenterOfMass.Velocity = COMV
        self.CenterOfMass.AngularMomentum = AngMom

        # now actually do the orientation
        AngMom_dir = AngMom/np.linalg.norm(AngMom)
        theta = np.arccos(np.dot(AngMom_dir, np.array([0, 0, 1])))
        k = np.cross(AngMom, np.array([0, 0, 1.]))
        k /= np.linalg.norm(k)

        self.CenterOfMass.theta = theta
        self.CenterOfMass.k = k
        
        for pt,Npart in enumerate(self.sn.NumPart_Total):
            if Npart == 0:
                continue
            
            part = getattr(self.sn, 'part'+str(pt))
            pos = part.Coordinates.value - COM
            vel = part.Velocities.value - COMV
    
            pos_rot = rodrigues_formula(k, pos.astype(np.float64), theta)
            vel_rot = rodrigues_formula(k, vel.astype(np.float64), theta)
        
            part.RotatedCoordinates = np.copy(pos_rot)
            part.RotatedVelocities = np.copy(vel_rot)
        
            part.rotpos = part.RotatedCoordinates.view()
            part.rotvel = part.RotatedVelocities.view()
                