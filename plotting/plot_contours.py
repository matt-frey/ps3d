from utils import *
from nc_reader import nc_reader
import matplotlib.pyplot as plt

import argparse
import os
from mpl_toolkits.axes_grid1 import ImageGrid
import colorcet as cc

parser = argparse.ArgumentParser(description='Create contour plot.')
parser.add_argument('--filename',
                    type=str,
                    help='output file')
parser.add_argument('--step',
                    type=int,
                    help='step number')
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
step = args.step
save_path = args.save_path
overwrite = args.overwrite
fignum = args.fignum

ncr = nc_reader()

ncr.open(fname)

t = ncr.get_all('t')

ix = 32 # - 3\pi/8
iys = 64  # -pi/4
iye = 145 # it is actually 144 due to the : operator, hence \pi/16
izs = 240 # 7\pi/16
ize = 257 # it is actually 256 due to the : operator, hence \pi/2

u = ncr.get_dataset(step=step, name='x_velocity')
v = ncr.get_dataset(step=step, name='y_velocity')
w = ncr.get_dataset(step=step, name='z_velocity')
zeta = ncr.get_dataset(step=step, name='z_vorticity')

xg, yg, zg = ncr.get_meshgrid()

ncr.close()


fig = plt.figure(figsize=(8, 5), dpi=200)
grid = ImageGrid(fig, 111,
                nrows_ncols=(2, 1),
                aspect=False,
                axes_pad=(0.1, 0.35),
                direction='row',
                share_all=True,
                cbar_location="right",
                cbar_mode='each',
                cbar_size="2%",
                 cbar_pad=0.05)

for i in range(2):
    grid[i].quiver(yg[ix, iys:iye, izs:ize],
                   zg[ix, iys:iye, izs:ize],
                   v[ix, iys:iye, izs:ize],
                   w[ix, iys:iye, izs:ize],
                   units='xy')

cf1 = grid[0].contourf(yg[ix, iys:iye, izs:ize], zg[ix, iys:iye, izs:ize],
                       zeta[ix, iys:iye, izs:ize], levels=50,
                       cmap=cc.cm['rainbow4_r'], zorder=-2)

cf2 = grid[1].contourf(yg[ix, iys:iye, izs:ize], zg[ix, iys:iye, izs:ize],
                       u[ix, iys:iye, izs:ize], levels=50,
                       cmap=cc.cm['rainbow4_r'], zorder=-2)

cbar1 = grid[0].cax.colorbar(cf1)
cbar1.formatter.set_powerlimits((0, 0))
cbar1.ax.yaxis.set_offset_position('left')

cbar2 = grid[1].cax.colorbar(cf2)
cbar2.formatter.set_powerlimits((0, 0))
cbar2.ax.yaxis.set_offset_position('left')


for i in range(2):
    grid[i].set_yticks(np.pi*np.array([7/16, 15/32, 1/2]),
                       [r'$7\pi/16$', r'$15\pi/32$', r'$\pi/2$'])

    grid[i].set_xticks(np.pi*np.array([1/16, 0, -1/16, -1/8, -3/16, -1/4]),
                       [r'$\pi/16$', r'$0$', r'$-\pi/16$', r'$-\pi/8$', r'$-3\pi/16$', r'$-\pi/4$'])

    grid[i].set_xlabel(r'$y$')
    grid[i].set_ylabel(r'$z$')


remove_xticks(grid[0])

add_annotation(grid[0], r'vertical vorticity,  $\zeta$', xy=(0.01, 1.07), fontsize=12)
add_annotation(grid[1], r'horizontal velocity, $u$', xy=(0.01, 1.07), fontsize=12)


add_annotation(grid[0], r'$x = -3\pi/8$', xy=(-0.12, 1.25), fontsize=12)
add_timestamp(grid[0], t[step], xy=(0.05, 1.25), fontsize=12)
    
save_figure(plt=plt, figpath=save_path, fignum=fignum, overwrite=overwrite)
plt.close()
