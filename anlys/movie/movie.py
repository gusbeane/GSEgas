import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import h5py as h5
from tqdm import tqdm
import sys

from matplotlib.animation import FuncAnimation

time_conv = 977.793 # converts time units to Myr

class animate_maker(object):
    def __init__(self, H, time, im, txt, vmin, vmax):
        self.H = H
        self.time = time
        
        self.im = im
        self.txt = txt

        self.vmin = vmin
        self.vmax = vmax

    def __call__(self, frame):
        self.im.set_data(self.H[frame].T)
        
        self.txt.set_text('t='+'{0:.2f}'.format(self.time[frame])+' Gyr')
        
        return (self.im, self.txt)

def make_movie(H, time, nres, vmin, vmax, fout,
               fps=16):
    
    # initialize fig and ax, remove boundary
    fig, ax = plt.subplots(1, figsize=(2, 2))
    fig.subplots_adjust(0, 0, 1, 1)
    ax.axis("off")

    # initialize im
    im = ax.imshow(np.full((nres, nres), vmin), origin='lower', 
                   norm=mpl.colors.LogNorm(vmin=vmin, vmax=vmax))

    # initialize time if needed
    txt = ax.text(0.6, 0.85, ' ', transform=ax.transAxes, fontsize=7, c='k')

    animate = animate_maker(H, time, im, txt, vmin, vmax)


    # initialize animator
    animation = FuncAnimation(fig, animate, tqdm(range(len(H))), interval=1000 / fps)
    animation.save(fout, dpi=nres)

if __name__ == '__main__':
    import h5py as h5

    t = h5.File('proj-SMUGGLE-point-lvl3.h5', mode='r')
    H = t['Hxy_g']
    time = np.arange(0, len(H)/200., 0.005)

    make_movie(H, time, H.shape[1], 1E-4, 1E-1, 'test.mp4')
