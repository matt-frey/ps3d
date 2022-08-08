from nc_reader import nc_reader
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
from utils import *
import argparse
import os

parser = argparse.ArgumentParser(description='Plot vorticity, velocity and helicty evolution.')
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

parser.add_argument('--fignums',
                    type=int,
                    nargs=2,
                    help='figure number')

args = parser.parse_args()
fname = args.filename
restart_file = args.restartfile
save_path = args.save_path
overwrite = args.overwrite
fignums = args.fignums

print()
print("\tFilename:    ", fname)
print("\tRestart file:", restart_file)
print("\tSave path:   ", save_path)
print("\tOverwrite:   ", overwrite)
print("\tFignums:     ", fignums)
print()

def fill_steps(ncr, j, lo, hi):

    t_all = ncr.get_all('t')

    for step in range(lo, hi):
        x_vor = ncr.get_dataset(step, 'x_vorticity')
        y_vor = ncr.get_dataset(step, 'y_vorticity')
        z_vor = ncr.get_dataset(step, 'z_vorticity')

        x_vel = ncr.get_dataset(step, 'x_velocity')
        y_vel = ncr.get_dataset(step, 'y_velocity')
        z_vel = ncr.get_dataset(step, 'z_velocity')

        H = ncr.get_dataset(step, 'helicity')
        he[j] = H.mean(axis=(0, 1, 2))

        nz, ny, nx = x_vor.shape

        t[j] = t_all[step]

        # add +1 to nz to include nz in the list
        izs = np.arange(0, nz+1, int(nz/4))
        for i, iz in enumerate(izs):
            u_rms[j, i]    = np.sqrt((x_vel[:, :, iz] ** 2).sum() / (nx * ny))
            v_rms[j, i]    = np.sqrt((y_vel[:, :, iz] ** 2).sum() / (nx * ny))
            w_rms[j, i]    = np.sqrt((z_vel[:, :, iz] ** 2).sum() / (nx * ny))
            xi_rms[j, i]   = np.sqrt((x_vor[:, :, iz] ** 2).sum() / (nx * ny))
            eta_rms[j, i]  = np.sqrt((y_vor[:, :, iz] ** 2).sum() / (nx * ny))
            zeta_rms[j, i] = np.sqrt((z_vor[:, :, iz] ** 2).sum() / (nx * ny))
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

t = np.zeros(n)
he = np.zeros(n)
u_rms = np.zeros((n, 5))
v_rms = np.zeros((n, 5))
w_rms = np.zeros((n, 5))

xi_rms = np.zeros((n, 5))
eta_rms = np.zeros((n, 5))
zeta_rms = np.zeros((n, 5))

if restart_file is None:
    fill_steps(ncr1, 0, 0, n)
else:
    fill_steps(ncr1, 0, 0, lo1)
    fill_steps(ncr2, lo1, 0, len(t2))
    fill_steps(ncr1, lo1+len(t2), hi1, len(t1))
    ncr2.close()

ncr1.close()

restart_label = ['restart', None]


#
# Helicity:
#
fig, ax = plt.subplots(1, 1, figsize=(8, 2), dpi=200)

ax.axvspan(xmin=t2[0], xmax=t2[-1], color='lightgrey', zorder=-1, label=restart_label[0])
ax.plot(t, he, label=r'$\mathcal{H}(t)$', color=colors[0])
ax.axhline(he[0], color='black', linestyle='dashdot', label=r'$\mathcal{H}(0)$')
ax.legend(loc='upper center', ncol=6, bbox_to_anchor=(0.5, 1.32))
ax.set_xlim([-1, 101])
ax.grid(zorder=-2)
ax.set_xlabel(r'time, $t$')

plt.tight_layout()
save_figure(plt=plt, figpath=save_path, fignum=fignums[0], overwrite=overwrite)
plt.close()


#
# Vorticity and velocity:
#
fig, axs = plt.subplots(2, 1, figsize=(8, 4), dpi=200, sharex=True, sharey=False)
grid = axs.flatten()

for k in range(2):
    grid[k].axvspan(xmin=t2[0], xmax=t2[-1], color='lightgrey', zorder=-1, label=restart_label[k])


# average lower and upper surface values:
u_rms[:, 0] = 0.5 * (u_rms[:, 0] + u_rms[:, 4])
v_rms[:, 0] = 0.5 * (v_rms[:, 0] + v_rms[:, 4])
w_rms[:, 0] = 0.5 * (w_rms[:, 0] + w_rms[:, 4])

grid[0].plot(t, u_rms[:, 0], label=r'$\langle \partial u_{\mathrm{rms}}\rangle$',
             color=colors[0], linestyle='dashed')
grid[0].plot(t, v_rms[:, 0], label=r'$\langle \partial v_{\mathrm{rms}}\rangle$',
             color=colors[1], linestyle='dashed')
grid[0].plot(t, w_rms[:, 0], label=r'$\langle\partial w_{\mathrm{rms}}\rangle$',
             color=colors[2], linestyle='dashed')

# middle
avg = (u_rms[:, 1] + u_rms[:, 2] + u_rms[:, 3]) / 3
grid[0].plot(t, avg, label=r'$\langle u_{\mathrm{rms}}\rangle$', color=colors[0],
             linestyle='solid')
avg = (v_rms[:, 1] + v_rms[:, 2] + v_rms[:, 3]) / 3
grid[0].plot(t, avg, label=r'$\langle v_{\mathrm{rms}}\rangle$', color=colors[1],
             linestyle='solid')
avg = (w_rms[:, 1] + w_rms[:, 2] + w_rms[:, 3]) / 3
grid[0].plot(t, avg, label=r'$\langle w_{\mathrm{rms}}\rangle$', color=colors[2],
             linestyle='solid')


# average lower and upper surface values:
xi_rms[:, 0] = 0.5 * (xi_rms[:, 0] + xi_rms[:, 4])
eta_rms[:, 0] = 0.5 * (eta_rms[:, 0] + eta_rms[:, 4])
zeta_rms[:, 0] = 0.5 * (zeta_rms[:, 0] + zeta_rms[:, 4])

grid[1].plot(t, xi_rms[:, 0], label=r'$\langle\partial\xi_{\mathrm{rms}}\rangle$',
             color=colors[0], linestyle='dashed')
grid[1].plot(t, eta_rms[:, 0], label=r'$\langle\partial\eta_{\mathrm{rms}}\rangle$',
             color=colors[1], linestyle='dashed')
grid[1].plot(t, zeta_rms[:, 0], label=r'$\langle\partial\zeta_{\mathrm{rms}}\rangle$',
             color=colors[2], linestyle='dashed')

# middle
avg = (xi_rms[:, 1] + xi_rms[:, 2] + xi_rms[:, 3]) / 3
grid[1].plot(t, avg, label=r'$\langle\xi_{\mathrm{rms}}\rangle$', color=colors[0],
             linestyle='solid')
avg = (eta_rms[:, 1] + eta_rms[:, 2] + eta_rms[:, 3]) / 3
grid[1].plot(t, avg, label=r'$\langle\eta_{\mathrm{rms}}\rangle$', color=colors[1],
             linestyle='solid')
avg = (zeta_rms[:, 1] + zeta_rms[:, 2] + zeta_rms[:, 3]) / 3
grid[1].plot(t, avg, label=r'$\langle\zeta_{\mathrm{rms}}\rangle$', color=colors[2],
             linestyle='solid')

for k in range(2):
    grid[k].legend(loc='upper center', ncol=6, bbox_to_anchor=(0.5, 1.32))
    grid[k].set_xlim([-1, 101])
    grid[k].grid(zorder=-2)
grid[1].set_xlabel(r'time, $t$')

# 3 August 2022
# https://stackoverflow.com/a/29988431
grid[0].tick_params(axis='x', which='both', length=0)

plt.tight_layout()
save_figure(plt=plt, figpath=save_path, fignum=fignums[1], overwrite=overwrite)
plt.close()
