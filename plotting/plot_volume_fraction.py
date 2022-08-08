from nc_reader import nc_reader
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
from utils import *
import argparse
import os

parser = argparse.ArgumentParser(description='Create cross section figures.')
parser.add_argument('--filename',
                    type=str,
                    help='output file')

parser.add_argument('--restartfile',
                    type=str,
                    help="another file restarted from the first one")

parser.add_argument('--save_path',
                    type=str,
                    help='where to save the figures',
                    default=os.getcwd())

parser.add_argument('--overwrite',
                    help='overwrite figures',
                    action='store_true')

parser.add_argument('--fignum',
                    type=int,
                    help='figure number')

args = parser.parse_args()
fname = args.filename
restart_file = args.restartfile
save_path = args.save_path
overwrite = args.overwrite
fignum = args.fignum

print()
print("\tFilename:    ", fname)
print("\tRestart file:", restart_file)
print("\tSave path:   ", save_path)
print("\tOverwrite:   ", overwrite)
print("\tFignum:      ", fignum)
print()

def fill_steps(ncr, j, lo, hi):

    t_all = ncr.get_all('t')

    for step in range(lo, hi):
        vor_mag = ncr.get_dataset(step, 'vorticity_magnitude')

        nx, ny, nz = vor_mag.shape

        #print("before:", nx, ny, nz)

        # remove periodic layers:
        vor_mag = vor_mag[0:nx-1, :, :]
        vor_mag = vor_mag[:, 0:ny-1, :]

        nx, ny, nz = vor_mag.shape

        #print("after:", nx, ny, nz)

        t[j] = t_all[step]

        rms = 0.5*(vor_mag[:,:,0]**2 + vor_mag[:,:,nz-1]**2).sum() + (vor_mag[:, :, 1:nz-1]**2).sum()

        rms = np.sqrt(rms / (nx * ny * (nz-1)))

        num = float((vor_mag > rms).sum())

        # volume cancels out:
        # ncell = nx * ny * (nz-1)
        # num * vcell / (total volume) = num * (pi**3 / ncell) / pi**3
        #                              = num / ncell
        #                              = num / (nx * ny * (nz-1))
        frac[j] = num / (nx * ny * (nz-1)) * 100.0 # in percent
        
        j = j + 1


ncr1 = nc_reader()
ncr2 = nc_reader()

ncr1.open(fname)

t1 = ncr1.get_all('t')
n = len(t1)

if not restart_file is None:
    ncr2.open(restart_file)
    t2 = ncr2.get_all('t')
    lo1 = find_nearest(t1, t2[0])
    hi1 = find_nearest(t1, t2[-1])
    n = len(t1[0:lo1]) + len(t2) + len(t1[hi1:])

fig, axs = plt.subplots(1, 1, figsize=(8, 3), dpi=200, sharex=True, sharey=False)
#grid = axs.flatten()


t = np.zeros(n)
frac = np.zeros(n)

if restart_file is None:
    fill_steps(ncr1, 0, 0, n)
else:
    fill_steps(ncr1, 0, 0, lo1)
    fill_steps(ncr2, lo1, 0, len(t2))
    fill_steps(ncr1, lo1+len(t2), hi1, len(t1))
    ncr2.close()

ncr1.close()

axs.axvspan(xmin=t2[0], xmax=t2[-1], color='lightgrey', zorder=-1, label='restart')


#axs.plot(t, he, label=r'$\mathcal{H}(t)$', color=colors[0])
#axs.axhline(he[0], color='black', linestyle='dashdot', label=r'$\mathcal{H}(0)$')
axs.plot(t, frac)  #, label=r'$\langle\zeta_{\mathrm{rms}}\rangle$', color=colors[2],
             #linestyle='solid')

#remove_xticks(grid[0])

#for k in range(3):
#    grid[k].legend(loc='upper center', ncol=6, bbox_to_anchor=(0.5, 1.32))
#    grid[k].set_xlim([-1, 101])
#    grid[k].grid(zorder=-2)
axs.set_xlabel(r'time, $t$')
axs.set_ylabel(r'volume fraction of $|\bm{\omega}|>|\bm{\omega}|_{\mathrm{rms}}$ ($\%$)')

# 3 August 2022
# https://stackoverflow.com/a/29988431
#grid[0].tick_params(axis='x', which='both', length=0)
#grid[1].tick_params(axis='x', which='both', length=0)

plt.tight_layout()
save_figure(plt=plt, figpath=save_path, fignum=fignum, overwrite=overwrite)
plt.close()
