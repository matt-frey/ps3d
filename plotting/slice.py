import numpy as np
import matplotlib.pyplot as plt
from tools.nc_reader import nc_reader
import matplotlib as mpl
import colorcet as cc
from mpl_toolkits.axes_grid1 import ImageGrid


def copy_periodic_layers(field):
    nz, ny, nx = field.shape
    field_copy = np.empty((nz, ny+1, nx+1))
    field_copy[:, 0:ny, 0:nx] = field.copy()
    field_copy[:, ny, :] = field_copy[:, 0, :]
    field_copy[:, :, nx] = field_copy[:, :, 0]
    return field_copy


def add_timestamp(plt, time, xy=(0.75, 1.05), fmt="%.3f"):
    # 29. Dec 2020
    # https://matplotlib.org/3.1.1/gallery/pyplots/annotate_transform.html#sphx-glr-gallery-pyplots-annotate-transform-py
    # https://stackoverflow.com/questions/7045729/automatically-position-text-box-in-matplotlib
    # https://matplotlib.org/3.1.0/gallery/recipes/placing_text_boxes.html
    bbox = dict(boxstyle="round", facecolor="wheat", edgecolor='none') #, alpha=0.5

    label = r"t = " + fmt % (time)

    plt.annotate(
        label, xy=xy, xycoords="axes fraction", bbox=bbox
    )

grid = '32'

mpl.rcParams['font.size'] = 10

ncreader = nc_reader()
ncreader.open('beltrami_' + grid + '_fields.nc')

z = ncreader.get_all('z')
t = ncreader.get_all('t')
origin = ncreader.get_box_origin()
extent = ncreader.get_box_extent()
ncells = ncreader.get_box_ncells()

xmin = origin[0]
xmax = origin[0] + extent[0]
ymin = origin[1]
ymax = origin[1] + extent[1]
zmin = origin[2]
zmax = origin[2] + extent[2]


steps = [60, 62, 64, 66, 68, 70]
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
ticks   = np.pi * np.array([-0.5, -0.25, 0.0, 0.25, 0.5])
ticklab = [r'$-\pi/2$', r'$-\pi/4$', r'$0$', r'$\pi/4$', r'$\pi/2$']
for i, step in enumerate(steps):
    x_vor = ncreader.get_dataset(step=step, name='x_vorticity')
    y_vor = ncreader.get_dataset(step=step, name='y_vorticity')
    z_vor = ncreader.get_dataset(step=step, name='z_vorticity')

    vor = np.sqrt(x_vor ** 2 + y_vor ** 2 + z_vor ** 2)
    vor = copy_periodic_layers(vor)

    print(t[step])

    ax = grid[i]
    im = ax.imshow(vor[iz, :, :], cmap=cc.cm['rainbow4'],
                   interpolation='bilinear',
                   origin='lower',
                   extent=[xmin, xmax, ymin, ymax])
    cbar = ax.cax.colorbar(im)

    if i < 3:
        ax.tick_params(
            axis='x',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            bottom=False,      # ticks along the bottom edge are off
            top=False,         # ticks along the top edge are off
            labelbottom=False) # labels along the bottom edge are off

    if i == 0 and i == 3:
        pass
    else:
        ax.tick_params(
            axis='y',
            which='both',
            right=False,
            left=False)

    ax.set_xticks(ticks, ticklab)
    ax.set_yticks(ticks, ticklab)
    ax.set_xlabel(r'$x$')
    ax.set_ylabel(r'$y$')

    add_timestamp(ax, t[step], xy=(0.03, 1.05), fmt="%.2f")

bbox = dict(boxstyle="round", facecolor="wheat", edgecolor='none')
grid[2].annotate(r'$z = -\pi/2$', xy=(0.6, 1.2), xycoords="axes fraction", bbox=bbox, fontsize=12)

plt.savefig('lower_surface_vorticity_magnitude.eps', format='eps')
plt.close()

#
# Pressure
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

for i, step in enumerate(steps):
    pres = ncreader.get_dataset(step=step, name='pressure')
    pres = copy_periodic_layers(pres)

    print(t[step])

    ax = grid[i]
    im = ax.imshow(pres[iz, :, :], cmap=cc.cm['rainbow4'],
                   interpolation='bilinear',
                   origin='lower',
                   extent=[xmin, xmax, ymin, ymax]) #,
                   #vmin=pmin,
                   #vmax=pmax)
    cbar = ax.cax.colorbar(im)

    if i < 3:
        ax.tick_params(
            axis='x',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            bottom=False,      # ticks along the bottom edge are off
            top=False,         # ticks along the top edge are off
            labelbottom=False) # labels along the bottom edge are off

    if i == 0 and i == 3:
        pass
    else:
        ax.tick_params(
            axis='y',
            which='both',
            right=False,
            left=False)

    ax.set_xticks(ticks, ticklab)
    ax.set_yticks(ticks, ticklab)
    ax.set_xlabel(r'$x$')
    ax.set_ylabel(r'$y$')

    add_timestamp(ax, t[step], xy=(0.03, 1.05), fmt="%.2f")

bbox = dict(boxstyle="round", facecolor="wheat", edgecolor='none')
grid[2].annotate(r'$z = -\pi/2$', xy=(0.6, 1.2), xycoords="axes fraction", bbox=bbox, fontsize=12)

plt.savefig('lower_surface_pressure.eps', format='eps')
plt.close()

ncreader.close()
