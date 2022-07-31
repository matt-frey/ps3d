import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from nc_reader import nc_reader

grid = 32
mpl.rcParams['font.size'] = 12

def add_timestamp(plt, time, xy=(0.75, 1.05), fmt="%.3f"):
    # 29. Dec 2020
    # https://matplotlib.org/3.1.1/gallery/pyplots/annotate_transform.html#sphx-glr-gallery-pyplots-annotate-transform-py
    # https://stackoverflow.com/questions/7045729/automatically-position-text-box-in-matplotlib
    # https://matplotlib.org/3.1.0/gallery/recipes/placing_text_boxes.html
    bbox = dict(boxstyle="round", facecolor="wheat", alpha=0.5)

    label = r"t = " + fmt % (time)

    plt.annotate(
        label, xy=xy, xycoords="axes fraction", bbox=bbox,
        fontsize=10
    )


def get_screenshot(grid, step):
    ncreader = nc_reader()
    ncreader.open('beltrami_' + str(grid) + '_fields.nc')
    t = ncreader.get_all('t')
    origin = ncreader.get_box_origin()
    extent = ncreader.get_box_extent()
    ncells = ncreader.get_box_ncells()
    x_vor = ncreader.get_dataset(step=step, name='x_vorticity')
    y_vor = ncreader.get_dataset(step=step, name='y_vorticity')
    z_vor = ncreader.get_dataset(step=step, name='z_vorticity')

    nz, ny, nx = x_vor.shape

    vor = np.sqrt(x_vor ** 2 + y_vor ** 2 + z_vor ** 2)

    ncreader.close()

    xmin = origin[0]
    xmax = origin[0] + extent[0]
    ymin = origin[1]
    ymax = origin[1] + extent[1]
    zmin = origin[2]
    zmax = origin[2] + extent[2]

    mlab.options.offscreen = True
    mlab.figure(size=(256, 256), bgcolor=(1,1,1), fgcolor=(0.,0.,0.))
    sf = mlab.pipeline.scalar_field(vor)
    obj = mlab.pipeline.volume(sf)
    mlax = mlab.axes(xlabel='x', x_axis_visibility=True,
                ylabel='y', y_axis_visibility=True,
                zlabel='z', z_axis_visibility=True,
                ranges=[xmin, xmax, ymin, ymax, zmin, zmax])

    mlax.title_text_property.bold = False
    mlax.title_text_property.italic = True
    mlax.title_text_property.font_family = 'arial'
    mlax.label_text_property.bold = False
    mlax.label_text_property.italic = True
    mlax.label_text_property.font_family = 'arial'
    #mlax.axes.font_factor = 2.0


    #colorbar = mlab.colorbar(object=obj,
                #label_fmt='%.2f',
                #orientation='vertical')

    ## 1 July 2022
    ## https://github.com/enthought/mayavi/issues/643
    #colorbar.scalar_bar.unconstrained_font_size = True
    #colorbar.label_text_property.font_family = 'arial'
    #colorbar.label_text_property.bold = False
    #colorbar.label_text_property.font_size=24
    #colorbar.scalar_bar_representation.position = (0.9, 0.05)
    mlab.view(azimuth=40.0, elevation=60, distance='auto', focalpoint='auto')
    return t[step], mlab.screenshot()

# 7 July 2022
# https://stackoverflow.com/questions/2566412/find-nearest-value-in-numpy-array
def find_nearest(t, tref):
    idx = (np.abs(tref - t)).argmin()
    print('Closest time to', tref, 'is', t[idx])
    return idx



def add_inset(ax, bounds, t, ke, step):

    tref, arr = get_screenshot(grid, step)

    idx = find_nearest(t=t, tref=tref)

    axins = ax.inset_axes(bounds=bounds, zorder=-1)
    ax.indicate_inset(bounds=[t[idx], ke[idx], 0.0, 0.0],
                    inset_ax=axins, zorder=1,
                    edgecolor='darkgrey', alpha=0.8)
    axins.imshow(arr)
    add_timestamp(axins, t[idx], xy=(0.05, 0.9), fmt="%.2f")

    # 7 July 2022
    # https://stackoverflow.com/questions/37039685/hide-tick-label-values-but-keep-axis-labels
    axins.set_xticks([])
    axins.set_yticks([])

    # 7 July 2022
    # https://stackoverflow.com/questions/1982770/matplotlib-changing-the-color-of-an-axis
    axins.spines['bottom'].set_color('darkgrey')
    axins.spines['top'].set_color('darkgrey')
    axins.spines['right'].set_color('darkgrey')
    axins.spines['left'].set_color('darkgrey')


mpl.rcParams.update({
    "figure.figsize": (9, 6),
    "figure.dpi": 200,
    "font.family": "serif",
    "font.size": 20,
    "text.usetex": True,
    'text.latex.preamble': "\n".join([
        r"\usepackage{amsmath}",
        r"\usepackage[utf8]{inputenc}",
        r"\usepackage[T1]{fontenc}",
#        r"\usepackage{siunitx}",
        ])
})


plt.figure(figsize=(13, 6.5), dpi=200)
ax = plt.gca()

t, ke, en = np.loadtxt('beltrami_' + str(grid) + '_ecomp.asc',
                       skiprows=1, unpack=True)

ncelli = 1.0 / grid ** 3

# calculate mean KE and mean EN
ke *= ncelli
ax.plot(t, ke)
ax.set_xlabel(r'time, $t$')
ax.set_ylabel(r'average kinetic energy, $\langle\mathcal{K}\rangle$')

w = 0.45
h = w
add_inset(ax=ax, bounds=[-0.1, 0.47, w, h], t=t, ke=ke, step=0)
add_inset(ax=ax, bounds=[ 0.17, 0.47, w, h], t=t, ke=ke, step=60)
add_inset(ax=ax, bounds=[ 0.17, 0.01, w, h], t=t, ke=ke, step=66)
add_inset(ax=ax, bounds=[ 0.6, 0.47, w, h], t=t, ke=ke, step=100)

plt.tight_layout()

#plt.show()

plt.savefig('figxx1.eps', format='eps') #, bbox_inches='tight', dpi=200)

plt.close()
