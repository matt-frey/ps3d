import numpy as np
import matplotlib.pyplot as plt
from nc_reader import nc_reader
import matplotlib.colors as mpl_colors
from mpl_toolkits.axes_grid1 import ImageGrid
from utils import *
import argparse
import os

parser = argparse.ArgumentParser(description='Create cross section figure of exact Beltrami.')

parser.add_argument('--plane',
                    type=str,
                    help="'xy', 'xz' or 'yz'")

parser.add_argument('--ngrid',
                    type=int,
                    help='number of grid cells')

parser.add_argument('--loc',
                    type=int,
                    default=0,
                    help='index location to extract plane')

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

parser.add_argument('--zlabel',
                    default=r'$z = -\pi/2$',
                    help='z-label for location')

args = parser.parse_args()
plane = args.plane
loc = args.loc
save_path = args.save_path
overwrite = args.overwrite
fignum = args.fignum
zlabel = args.zlabel
ngrid = args.ngrid

print()
print("\tPlane:     ", plane)
print("\tLocation:  ", loc)
print("\tSave path: ", save_path)
print("\tOverwrite: ", overwrite)
print("\tFignum:   ", fignum)
print()

origin = -0.5 * np.pi * np.ones(3)
extent = np.pi * np.ones(3)

dx = extent / ngrid

k = 2.0
l = 2.0
m = 1.0

alpha = np.sqrt(k ** 2 + l ** 2 + m ** 2)

nz = ngrid
ny = ngrid
nx = ngrid


x = np.linspace(origin[0], origin[0]+extent[0], nx)
y = np.linspace(origin[1], origin[1]+extent[1], ny)
z = np.linspace(origin[2], origin[2]+extent[2], nz)

xg, yg, zg = np.meshgrid(x, y, z, indexing='ij')

x_vel = np.zeros((nz, ny, nx))
y_vel = np.zeros((nz, ny, nx))
z_vel = np.zeros((nz, ny, nx))
x_vor =	np.zeros((nz, ny, nx))
y_vor =	np.zeros((nz, ny, nx))
z_vor =	np.zeros((nz, ny, nx))

for ix in range(ngrid):
    for iy in range(ngrid):
        for iz in range(ngrid):
            x = origin[0] + ix * dx[0]
            y = origin[1] + iy * dx[1]
            z = origin[2] + iz * dx[2]
            
            u = 1.0 / (k**2 + l**2) * (k*m*np.sin(m*z) - l*alpha*np.cos(m*z)) * np.sin(k*x + l*y)
            v = 1.0 / (k**2 + l**2) * (l*m*np.sin(m*z) + k*alpha*np.cos(m*z)) * np.sin(k*x + l*y)
            w = np.cos(m*z) * np.cos(k*x + l*y)
        
            
            xi = alpha * u
            eta = alpha * v
            zeta = alpha * w
                
            x_vor[ix, iy, iz] = xi
            y_vor[ix, iy, iz] = eta
            z_vor[ix, iy, iz] = zeta
            x_vel[ix, iy, iz] = u
            y_vel[ix, iy, iz] = v
            z_vel[ix, iy, iz] = w

print("Initialized fields.")
            
fig = plt.figure(figsize=(8, 5), dpi=200)
grid = ImageGrid(fig, 111,
                nrows_ncols=(1, 3),
                aspect=True,
                axes_pad=(0.55, 0.3),
                direction='row',
                share_all=True,
                cbar_location="right",
                cbar_mode='each',
                cbar_size="4%",
                cbar_pad=0.05)

fields = [x_vor, y_vor, z_vor]
labels = [r'$\xi$', r'$\eta$', r'$\zeta$']

n_levels = 20

for i, field in enumerate(fields):
    ax = grid[i]
    # plane = 'xy':
    imin = origin[0]
    imax = origin[0] + extent[0]
    jmin = origin[1]
    jmax = origin[1] + extent[1]


    if plane == 'xy':
        pl = field[:, :, loc]
        xlab = r'$x$'
        ylab = r'$y$'
        xx = xg[:, :, loc]
        yy = yg[:, :, loc]
    elif plane == 'xz':
        jmin = origin[2]
        jmax = origin[2] + extent[2]
        pl = field[:, loc, :]
        xlab = r'$x$'
        ylab = r'$z$'
        xx = xg[:, loc, :]
        yy = zg[:, loc, :]
    elif plane == 'yz':
        imin = jmin
        imax = jmax
        jmin = origin[2]
        jmax = origin[2] + extent[2]
        pl = field[loc, :, :]
        xlab = r'$y$'
        ylab = r'$z$'
        xx = yg[loc, :, :]
        yy = zg[loc, :, :]

#    ax.quiver(xx, yy,
#              v[ix, iys:iye, izs:ize],
#              w[ix, iys:iye, izs:ize],
#              units='xy')

    cf1 = ax.contourf(xx, yy, pl, levels=n_levels,
                      cmap=cc.cm['coolwarm'], zorder=-2)

    cbar = ax.cax.colorbar(cf1)
    cbar.formatter.set_powerlimits((0, 0))
    cbar.ax.yaxis.set_offset_position('left')

    ticks   = np.pi * np.array([-0.5, -0.25, 0.0, 0.25, 0.5])
    ticklab = [r'$-\pi/2$', r'$-\pi/4$', r'$0$', r'$\pi/4$', r'$\pi/2$']
    ax.set_xticks(ticks, ticklab)
    ax.set_yticks(ticks, ticklab)
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)

    if i == 0 or i == 3:
        pass
    else:
        remove_yticks(ax)

    add_annotation(ax, labels[i], xy=(0.03, 1.06))

if not zlabel is None:
    add_annotation(grid[0], zlabel, xy=(-0.5, 1.25), fontsize=12)

save_figure(plt=plt, figpath=save_path, fignum=fignum, overwrite=overwrite)
plt.close()
