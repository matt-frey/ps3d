import numpy as np
import matplotlib.pyplot as plt
from utils import *
import argparse
import os
from numpy.polynomial import Polynomial as poly
import matplotlib.gridspec as gridspec

parser = argparse.ArgumentParser(
    description='Plot kinetic energy and enstrophy for different grid resolutions.')

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
                    default=3,
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

labels = [r'$32^3$',
          r'$64^3$',
          r'$128^3$',
          r'$256^3$'
]

grids = np.array([32, 64, 128, 256])

#fig, axs = plt.subplots(2, 1, figsize=(4, 4), dpi=200, sharex=True)

fig = plt.figure(figsize=(9, 4.5), dpi=200, tight_layout=True)
gs = gridspec.GridSpec(2, 2)

i = 0

ax1 = fig.add_subplot(gs[1, 0])
ax0 = fig.add_subplot(gs[0, 0]) #, sharex=ax1)
ax2 = fig.add_subplot(gs[:, 1])

for grid in grids:
    t, ke, en = np.loadtxt(os.path.join(fpath, 'beltrami_' + str(grid) + '_ecomp.asc'),
                           skiprows=1, unpack=True)

    voli = 1.0 / np.pi ** 3
    ncelli = 1.0 / grid ** 3

    # calculate mean KE and mean EN
    ke *= voli * ncelli
    en *= voli * ncelli

    #print("initial <KE>", ke[0])
    #print("initial <EN>", en[0])

    print("Initial KE of grid", grid, "is", ke[0])
    print("Final KE of grid", grid, "is", ke[-1])
    print("Remaining KE in percent:", ke[-1]/ke[0] * 100.0)

    label = labels[i]
    i = i + 1

    ax0.plot(t, ke / ke[0], label=label)
    ax1.plot(t, en / en[0], label=label)

ax1.set_xlabel(r'time, $t$')
ax1.set_ylabel(r'enstrophy, $\Upsilon(t)/\Upsilon(0)$')

ax0.tick_params(axis='x', which='both', length=0)
ax0.set_xticklabels([])

ax0.grid(zorder=-1)
ax1.grid(zorder=-1)

ax0.set_xlim([50, 101])
ax1.set_xlim([50, 101])

ax0.set_ylabel(r'kinetic energy, $\mathcal{K}(t)/\mathcal{K}(0)$')

ax0.legend(loc='upper center', ncol=5, bbox_to_anchor=(0.5, 1.3))



maxen = np.zeros(len(grids))
vmax = np.zeros(len(grids))

#fig, axs = plt.subplots(1, 1, figsize=(4, 4), dpi=200, sharex=False)

for i, grid in enumerate(grids):
    _, _, en = np.loadtxt(os.path.join(fpath, 'beltrami_' + str(grid) + '_ecomp.asc'),
                          skiprows=1, unpack=True)

    # t max rms char <xi> <eta> <zeta>
    _, vormax, _, _, _, _, _ = np.loadtxt(
        os.path.join(fpath, 'beltrami_' + str(grid) + '_vorticity.asc'),
        skiprows=1, unpack=True)

    ncelli = 1.0 / grid ** 3
    voli = 1.0 / np.pi ** 3

    # calculate mean KE and mean EN
    en *= ncelli * voli

    maxen[i] = en.max()
    vmax[i] = vormax.max()

# ignore nz = 32
log10_maxen = np.log10(maxen[1:])
log10_vormax = np.log10(vmax[1:])
log10_nz = np.log10(grids[1:])


# linear fit: log10(maxen) = m * log10(nz) + q
p_fitted = poly.fit(x=log10_nz, y=log10_maxen, deg=1) #, domain=[6, 8], window=[6, 8])
p_fitted = p_fitted.convert() # to unscale, i.e. makes domain == window --> same result as np.polyfit

np.polynomial.set_default_printstyle('ascii')
print("Fitted polynomial:", p_fitted)
q = p_fitted.coef[0]
m = p_fitted.coef[1]

#print("polyfit:", np.polyfit(x=log10_nz, y=log10_maxen, deg=1))
#print(p_fitted.coef)

ax2.plot(grids, maxen, marker='o', markersize=4, label=r'$\Upsilon$')
ax2.plot(grids, vmax, marker='o', markersize=4, label=r'$|\bm{\omega}|_{\max}$')
ax2.plot(grids[1:], 10 ** (m * log10_nz + q), linestyle='dashed', color='black',
         label=r'$\log_{10}\Upsilon\propto' + str(round(m, 3)) + '\log_{10}n_{z}$')


# linear fit: log10(vmax) = m * log10(nz) + q
p_fitted = poly.fit(x=log10_nz, y=log10_vormax, deg=1)
p_fitted = p_fitted.convert() # to unscale, i.e. makes domain == window --> same result as np.polyfit
np.polynomial.set_default_printstyle('ascii')
print("Fitted polynomial:", p_fitted)
q = p_fitted.coef[0]
m = p_fitted.coef[1]
ax2.plot(grids[1:], 10 ** (m * log10_nz + q), linestyle='dashed', color='black')


ax2.set_xscale('log', base=2)
ax2.set_yscale('log', base=10)

ax2.grid(zorder=-1)

ax2.set_xlabel(r'grid resolution ($n_x = n_y = n_z$)')
ax2.set_ylabel(r'peak values')
ax2.legend(loc='upper center', ncol=5, bbox_to_anchor=(0.5, 1.15))

plt.tight_layout()

save_figure(plt=plt, figpath=save_path, fignum=fignum, overwrite=overwrite)
plt.close()
