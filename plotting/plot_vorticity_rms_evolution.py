from tools.nc_reader import nc_reader
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
save_path = args.save_path
overwrite = args.overwrite
fignum = args.fignum

print()
print("\tFilename:  ", fname)
print("\tSave path: ", save_path)
print("\tOverwrite: ", overwrite)
print("\tFignum:    ", fignum)
print()

mpl.rcParams['font.size'] = 10

fig, axs = plt.subplots(2, 1, figsize=(8, 5), dpi=400, sharex=True, sharey=False)
grid = axs.flatten()

ncreader = nc_reader()
ncreader.open(fname)

t = ncreader.get_all('t')

n = len(t)

u_rms = np.zeros((n, 5))
v_rms = np.zeros((n, 5))
w_rms = np.zeros((n, 5))

xi_rms = np.zeros((n, 5))
eta_rms = np.zeros((n, 5))
zeta_rms = np.zeros((n, 5))

for step in range(n):

    x_vor = ncreader.get_dataset(step, 'x_vorticity')
    y_vor = ncreader.get_dataset(step, 'y_vorticity')
    z_vor = ncreader.get_dataset(step, 'z_vorticity')

    x_vel = ncreader.get_dataset(step, 'x_velocity')
    y_vel = ncreader.get_dataset(step, 'y_velocity')
    z_vel = ncreader.get_dataset(step, 'z_velocity')

    nz, ny, nx = x_vor.shape

    # add +1 to nz to include nz in the list
    izs = np.arange(0, nz+1, int(nz/4))
    for i, iz in enumerate(izs):
        u_rms[step, i]    = np.sqrt((x_vel[iz, :, :] ** 2).sum() / (nx * ny))
        v_rms[step, i]    = np.sqrt((y_vel[iz, :, :] ** 2).sum() / (nx * ny))
        w_rms[step, i]    = np.sqrt((z_vel[iz, :, :] ** 2).sum() / (nx * ny))
        xi_rms[step, i]   = np.sqrt((x_vor[iz, :, :] ** 2).sum() / (nx * ny))
        eta_rms[step, i]  = np.sqrt((y_vor[iz, :, :] ** 2).sum() / (nx * ny))
        zeta_rms[step, i] = np.sqrt((z_vor[iz, :, :] ** 2).sum() / (nx * ny))

ncreader.close()

# bottom
grid[0].plot(t, u_rms[:, 0], label=r'$u_{\mathrm{rms}}(z = -\pi/2)$',
             color='blue', linestyle='dotted')
grid[0].plot(t, v_rms[:, 0], label=r'$v_{\mathrm{rms}}(z = -\pi/2)$',
             color='red', linestyle='dotted')
grid[0].plot(t, w_rms[:, 0], label=r'$w_{\mathrm{rms}}(z = -\pi/2)$',
             color='green', linestyle='dotted')

# top
grid[0].plot(t, u_rms[:, 4], label=r'$u_{\mathrm{rms}}(z = \pi/2)$',
             color='blue', linestyle='dashed')
grid[0].plot(t, v_rms[:, 4], label=r'$v_{\mathrm{rms}}(z = \pi/2)$',
             color='red', linestyle='dashed')
grid[0].plot(t, w_rms[:, 4], label=r'$w_{\mathrm{rms}}(z = \pi/2)$',
             color='green', linestyle='dashed')

# middle
avg = (u_rms[:, 1] + u_rms[:, 2] + u_rms[:, 3]) / 3
grid[0].plot(t, avg, label=r'$\langle u_{\mathrm{rms}}\rangle$', color='orange',
             linestyle='solid')
avg = (v_rms[:, 1] + v_rms[:, 2] + v_rms[:, 3]) / 3
grid[0].plot(t, avg, label=r'$\langle v_{\mathrm{rms}}\rangle$', color='blue',
             linestyle='solid')
avg = (w_rms[:, 1] + w_rms[:, 2] + w_rms[:, 3]) / 3
grid[0].plot(t, avg, label=r'$\langle w_{\mathrm{rms}}\rangle$', color='red',
             linestyle='solid')


# bottom
grid[1].plot(t, xi_rms[:, 0], label=r'$\xi_{\mathrm{rms}}(z = -\pi/2)$',
             color='blue', linestyle='dotted')
grid[1].plot(t, eta_rms[:, 0], label=r'$\eta_{\mathrm{rms}}(z = -\pi/2)$',
             color='red', linestyle='dotted')
grid[1].plot(t, zeta_rms[:, 0], label=r'$\zeta_{\mathrm{rms}}(z = -\pi/2)$',
             color='green', linestyle='dotted')

# top
grid[1].plot(t, xi_rms[:, 4], label=r'$\xi_{\mathrm{rms}}(z = \pi/2)$',
            color='blue', linestyle='dashed')
grid[1].plot(t, eta_rms[:, 4], label=r'$\eta_{\mathrm{rms}}(z = \pi/2)$',
             color='red', linestyle='dashed')
grid[1].plot(t, zeta_rms[:, 4], label=r'$\zeta_{\mathrm{rms}}(z = \pi/2)$',
             color='green', linestyle='dashed')

# middle
avg = (xi_rms[:, 1] + xi_rms[:, 2] + xi_rms[:, 3]) / 3
grid[1].plot(t, avg, label=r'$\langle\xi_{\mathrm{rms}}\rangle$', color='orange',
             linestyle='solid')
avg = (eta_rms[:, 1] + eta_rms[:, 2] + eta_rms[:, 3]) / 3
grid[1].plot(t, avg, label=r'$\langle\eta_{\mathrm{rms}}\rangle$', color='blue',
             linestyle='solid')
avg = (zeta_rms[:, 1] + zeta_rms[:, 2] + zeta_rms[:, 3]) / 3
grid[1].plot(t, avg, label=r'$\langle\zeta_{\mathrm{rms}}\rangle$', color='red',
             linestyle='solid')

#remove_xticks(grid[0])

for k in range(2):
    grid[k].legend(loc='right', ncol=1, bbox_to_anchor=(1.28, 0.5))
    grid[k].set_xlim([-1, 101])
grid[1].set_xlabel(r'time, $t$')


plt.tight_layout()
save_figure(plt=plt, figpath=save_path, fignum=fignum, overwrite=overwrite)
plt.close()
