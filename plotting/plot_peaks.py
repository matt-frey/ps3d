import numpy as np
import matplotlib.pyplot as plt
from utils import *
import argparse
import os

parser = argparse.ArgumentParser(
    description='Plot kinetic energy and enstrophy for different hyperdiffusion prefactors.')

parser.add_argument('--filepath',
                    type=str,
                    help='path to files')

parser.add_argument('--save_path',
                    type=str,
                    help='where to save the figures',
                    default=os.getcwd())

parser.add_argument('--overwrite',
                    help='overwrite figures',
                    action='store_true')

parser.add_argument('--fignum',
                    type=int,
                    default=2,
                    help='figure number')

args = parser.parse_args()
fpath = args.filepath
save_path = args.save_path
overwrite = args.overwrite
fignum = args.fignum

print()
print("\tFilepath:  ", fpath)
print("\tSave path: ", save_path)
print("\tOverwrite: ", overwrite)
print("\tFignum:    ", fignum)
print()


grids = [32, 64, 128, 256]

maxen = np.zeros(len(grids))
vmax = np.zeros(len(grids))
vrms = np.zeros(len(grids))

mpl.rcParams['font.size'] = 11

fig, axs = plt.subplots(1, 1, figsize=(7, 2.5), dpi=200, sharex=True)

for i, grid in enumerate(grids):
    _, _, en = np.loadtxt(os.path.join(fpath, 'beltrami_' + str(grid) + '_ecomp.asc'),
                          skiprows=1, unpack=True)

    # t max rms char <xi> <eta> <zeta>
    _, vormax, vorrms, _, _, _, _ = np.loadtxt(os.path.join('beltrami_' + str(grid) + '_vorticity.asc'),
                                               skiprows=1, unpack=True)

    ncelli = 1.0 / grid ** 3

    # calculate mean KE and mean EN
    en *= ncelli

    maxen[i] = en.max()
    vmax[i] = vormax.max()
    vrms[i] = vorrms.max()

axs.plot(grids, maxen, marker='o', markersize=4, label=r'$\langle\Upsilon\rangle$')
axs.plot(grids, vmax, marker='o', markersize=4, label=r'$|\bm{\omega}|_{\max}$')
axs.plot(grids, vrms, marker='o', markersize=4, label=r'$|\bm{\omega}|_{\mathrm{rms}}$')
axs.set_xscale('log', base=2)
axs.set_yscale('log', base=10)

axs.grid(zorder=-1)

axs.set_xlabel(r'grid resolution')
axs.set_ylabel(r'peak values')
axs.legend(loc='upper center', ncol=5, bbox_to_anchor=(0.5, 1.35))

plt.tight_layout()

save_figure(plt=plt, figpath=save_path, fignum=fignum, overwrite=overwrite)
plt.close()
