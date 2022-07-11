import numpy as np
import matplotlib.pyplot as plt
from tools.nc_reader import nc_reader
import matplotlib as mpl
import colorcet as cc
from mpl_toolkits.axes_grid1 import ImageGrid

def add_timestamp(plt, time, xy=(0.75, 1.05), fmt="%.3f"):
    # 29. Dec 2020
    # https://matplotlib.org/3.1.1/gallery/pyplots/annotate_transform.html#sphx-glr-gallery-pyplots-annotate-transform-py
    # https://stackoverflow.com/questions/7045729/automatically-position-text-box-in-matplotlib
    # https://matplotlib.org/3.1.0/gallery/recipes/placing_text_boxes.html
    bbox = dict(boxstyle="round", facecolor="wheat", edgecolor='none') #, alpha=0.5,

    label = r"t = " + fmt % (time)

    plt.annotate(
        label, xy=xy, xycoords="axes fraction", bbox=bbox
    )

grid = '32'
step = 10

mpl.rcParams['font.size'] = 8

ncreader = nc_reader()
ncreader.open('beltrami_' + grid + '_fields.nc')

z = ncreader.get_all('z')
t = ncreader.get_all('t')
origin = ncreader.get_box_origin()
extent = ncreader.get_box_extent()
ncells = ncreader.get_box_ncells()

colors = ['blue', 'orange', 'green']

vcell = np.prod(extent / ncells)

steps = [0, 55, 60, 62, 64, 66, 68, 70, 80, 100]

n = len(z)
fig = plt.figure(figsize=(len(steps)*1.0, 2*3), dpi=400)


grid = ImageGrid(fig, 111,
                 nrows_ncols=(2, 5),
                 aspect=True,
                 axes_pad=(0.1, 0.3),
                 direction='row',
                 share_all=True,
                 cbar_location="bottom",
                 cbar_mode=None,
                 cbar_size="4%",
                 cbar_pad=0.0)
yticks   = np.pi * np.array([-0.5, -0.25, 0.0, 0.25, 0.5])
yticklab = [r'$-\pi/2$', r'$-\pi/4$', r'$0$', r'$\pi/4$', r'$\pi/2$']

for i, step in enumerate(steps):
    xibar = np.zeros(n)
    etabar = np.zeros(n)
    zetabar = np.zeros(n)
    #ubar = np.zeros(n)
    #vbar = np.zeros(n)
    #wbar = np.zeros(n)

    xi = ncreader.get_dataset(step=step, name='x_vorticity')
    eta = ncreader.get_dataset(step=step, name='y_vorticity')
    zeta = ncreader.get_dataset(step=step, name='z_vorticity')

    #u = ncreader.get_dataset(step=step, name='x_velocity')
    #v = ncreader.get_dataset(step=step, name='y_velocity')
    #w = ncreader.get_dataset(step=step, name='z_velocity')


    xibar = xi.mean(axis=(1, 2))
    etabar = eta.mean(axis=(1, 2))
    zetabar = zeta.mean(axis=(1, 2))

    #ubar = u.mean(axis=(1, 2))
    #vbar = v.mean(axis=(1, 2))
    #wbar = w.mean(axis=(1, 2))

    ax = grid[i]

    label1 = None
    label2 = None
    label3 = None
    if i == 2:
        label1 = r'$\langle\mathcal{\xi}\rangle$'
        label2 = r'$\langle\eta\rangle$'
        label3 = r'$\langle\zeta\rangle$'


    ax.plot(xibar, z, color='blue', marker='o', markersize=3, linewidth=0.75, label=label1)
    ax.plot(etabar, z, color='red', marker='x', markersize=3, linewidth=0.75, label=label2)
    ax.plot(zetabar, z, color='green', marker='+', markersize=3, linewidth=0.75, label=label3)
    ax.set_aspect(1)
    ax.set_yticks(ticks=yticks, labels=yticklab)
    #ax.set_xlabel(r'magn.')
    ax.set_xticks([-1, 0, 1])
    ax.set_xlim([-1.1, 1.1])

    add_timestamp(ax, t[step], xy=(0.03, 1.05), fmt="%.2f")

grid[0].set_ylabel(r'$z$')
grid[5].set_ylabel(r'$z$')

grid[2].legend(loc='upper center', ncol=3, bbox_to_anchor=(0.5, 1.25))

plt.tight_layout()
#plt.show()
plt.savefig('vor_profile.eps', dpi=400)

ncreader.close()
