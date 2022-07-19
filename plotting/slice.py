import numpy as np
import matplotlib.pyplot as plt
from nc_reader import nc_reader
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import ImageGrid
from utils import *

grid = '32'

ncreader = nc_reader()
ncreader.open('../examples/beltrami_' + grid + '_fields.nc')

t = ncreader.get_all('t')

steps = [0, 0, 0, 0, 0, 0]
iz = 0


#
# Vorticity magnitude
#
fig = plt.figure(figsize=(8, 5), dpi=200)
grid = ImageGrid(fig, 111,
                 nrows_ncols=(2, 3),
                 aspect=True,
                 axes_pad=(0.4, 0.3),
                 direction='row',
                 share_all=True,
                 cbar_location="right",
                 cbar_mode='each',
                 cbar_size="4%",
                 cbar_pad=0.1)

for i, step in enumerate(steps):
    x_vor = ncreader.get_dataset(step=step, name='x_vorticity')
    y_vor = ncreader.get_dataset(step=step, name='y_vorticity')
    z_vor = ncreader.get_dataset(step=step, name='z_vorticity')

    vor = np.sqrt(x_vor ** 2 + y_vor ** 2 + z_vor ** 2)

    axis = ncreader.get_axis('z')

    vor[:, :, :] = 0
    vor[:, 20, 0] = axis

    print(t[step])

    ax = grid[i]
    im, cbar = make_imshow(ax=ax,
                           plane='xy',
                           loc=iz,
                           fdata=vor,
                           ncr=ncreader,
                           cmap='rainbow4',
                           colorbar=True)

    if i < 3:
        remove_xticks(ax)

    if i == 0 and i == 3:
        pass
    else:
        remove_yticks(ax)

    add_timestamp(ax, t[step], xy=(0.03, 1.06), fmt="%.2f")

add_annotation(grid[2], r'$z = -\pi/2$', xy=(0.6, 1.2), fontsize=12)

plt.savefig('lower_surface_vorticity_magnitude.eps', format='eps')
plt.close()

#
# Pressure
#
fig = plt.figure(figsize=(8, 5), dpi=200)
grid = ImageGrid(fig, 111,
                 nrows_ncols=(1, 2),
                 aspect=True,
                 axes_pad=(0.4, 0.3),
                 direction='row',
                 share_all=True,
                 cbar_location="right",
                 cbar_mode='each',
                 cbar_size="4%",
                 cbar_pad=0.1)
ticks   = np.pi * np.array([-0.5, -0.25, 0.0, 0.25, 0.5])
ticklab = [r'$-\pi/2$', r'$-\pi/4$', r'$0$', r'$\pi/4$', r'$\pi/2$']

## find colorbar bounds
#pmin = 1000000
#pmax = -pmin
#for i, step in enumerate(steps):
    #pres = ncreader.get_dataset(step=step, name='pressure')
    #pres = copy_periodic_layers(pres)
    #pmin = min(pmin, pres[iz, :, :].min())
    #pmax = max(pmax, pres[iz, :, :].max())

for i, step in enumerate([0, 0]):
    pres = ncreader.get_dataset(step=step, name='pressure')

    X, Y, Z = ncreader.get_meshgrid()


    a = 2 * X + 2 * Y
    p = 0.25 * (0.125 * np.cos(2.0 * a) - np.cos(-np.pi))

    print (i)
    if i > 0:
        print ("HI")
        pres[:, :, iz] = p[:, :, iz]

    print(t[step])

    ax = grid[i]

    im, cbar = make_imshow(ax=ax,
                           plane='xy',
                           loc=iz,
                           fdata=pres,
                           ncr=ncreader,
                           cmap='rainbow4',
                           colorbar=True)
    if i < -1:
        remove_xticks(ax)

    if i == 0 and i == 3:
        pass
    else:
        remove_yticks(ax)

    ax.set_xticks(ticks, ticklab)
    ax.set_yticks(ticks, ticklab)
    ax.set_xlabel(r'$x$')
    ax.set_ylabel(r'$y$')

    #add_timestamp(ax, t[step], xy=(0.03, 1.05), fmt="%.2f")

add_annotation(grid[1], r'$z = -\pi/2$', xy=(0.6, 1.2), fontsize=12)

plt.savefig('lower_surface_pressure.eps', format='eps')
plt.close()

ncreader.close()
