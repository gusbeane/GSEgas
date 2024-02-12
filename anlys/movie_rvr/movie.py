import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import h5py as h5
from tqdm import tqdm
import sys

from matplotlib.animation import FuncAnimation

time_conv = 977.793 # converts time units to Myr

class animate_maker(object):
    def __init__(self, H, time, im, txt, norm):
        self.H = H
        self.time = time
        
        self.im = im
        self.txt = txt

        self.vmin = norm.vmin
        self.vmax = norm.vmax

    def __call__(self, frame):
        H_ = np.copy(self.H[frame])
        H_[H_ < self.vmin] = self.vmin
        #self.H[frame][self.H[frame] < self.vmin] = self.vmin

        self.im.set_data(H_.T)
        
        self.txt.set_text('t='+'{0:.2f}'.format(self.time[frame])+' Gyr')
        
        return (self.im, self.txt)

def make_movie(H, time, nresx, nresy, norm, fout, cmap, extent,
               fps=16):
    
    # initialize fig and ax, remove boundary
    fig, ax = plt.subplots(1, figsize=(4, 4))
    # fig.subplots_adjust(0, 0, 1, 1)
    # ax.axis("off")

    # initialize im
    im = ax.imshow(np.full((nresx, nresy), norm.vmin), origin='lower', 
                   norm=norm, cmap=cmap, extent=extent)

    ax.set(xlabel='r [kpc]', ylabel='v_r [km/s]')
    
    # initialize time if needed
    txt = ax.text(0.6, 0.85, ' ', transform=ax.transAxes, fontsize=12, c='k')
    
    fig.tight_layout()

    animate = animate_maker(H, time, im, txt, norm)


    # initialize animator
    animation = FuncAnimation(fig, animate, tqdm(range(len(H))), interval=1000 / fps)
    animation.save(fout, dpi=nresx)

if __name__ == '__main__':
    import h5py as h5

    t = h5.File('proj-SMUGGLE-point-lvl3.h5', mode='r')
    H = t['Hxy_g']
    time = np.arange(0, len(H)/200., 0.005)

    make_movie(H, time, H.shape[1], 1E-4, 1E-1, 'test.mp4')
