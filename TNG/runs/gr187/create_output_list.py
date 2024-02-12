import numpy as np

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

def get_scale_factor(time, Omega0=0.3089, OmegaLambda=0.6911, HubbleParam=0.6774):
    HUBBLE = 3.2407789e-18
    SEC_PER_MEGAYEAR = 3.15576e13

    ans = time * (HUBBLE * HubbleParam) * (SEC_PER_MEGAYEAR * 1000)

    fac = 2. / (3. * np.sqrt(OmegaLambda))
    ans /= fac

    ans = np.sinh(ans)
    ans = ans**2 * Omega0/OmegaLambda
    ans = ans**(1./3.)

    return ans

a0 = 1.0
t0 = get_time(a0)

snap_time = [a0]
snap_full = [1]

dt = 0.01
Nfull = 20
N = Nfull - 1
t = t0
while t > t0 - 12:
    t -= dt
    snap_time.append(get_scale_factor(t))
    if N == 0:
        snap_full.append(1)
        N = Nfull - 1
    else:
        snap_full.append(3)
        N -= 1

print('low time res starts at t=', t, 'z=', 1./get_scale_factor(t)-1)
        
dt_early = 0.1
Nfull = 5
N = 0

ai = 0.0078125
ti = get_time(ai)

t -= dt_early

while t > ti:
    snap_time.append(get_scale_factor(t))
    if N == 0:
        snap_full.append(1)
        N = Nfull - 1
    else:
        snap_full.append(3)
        N -= 1
    t -= dt_early

snap_time.append(ai)
snap_full.append(1)

out_arr = np.column_stack((snap_time, snap_full))
np.savetxt('output_list.txt', out_arr, fmt='%.16f %d')

print('tot number of snap=', len(out_arr))
print('tot number of full snaps=', np.sum(out_arr[:,1]==1))
print('tot number of mini snaps=', np.sum(out_arr[:,1]==3))

