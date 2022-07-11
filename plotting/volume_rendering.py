# mayavi2 -x volume_rendering.py -o
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
    bbox = dict(boxstyle="round", facecolor="wheat", alpha=0.5)

    label = r"t = " + fmt % (time)

    plt.annotate(
        label, xy=xy, xycoords="axes fraction", bbox=bbox
    )

grid = '32'
step = 10

mpl.rcParams['font.size'] = 18

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

colors = ['blue', 'orange', 'green']

vcell = np.prod(extent / ncells)

#pres = ncreader.get_dataset(step=step, name='pressure')
x_vor = ncreader.get_dataset(step=step, name='x_vorticity')
y_vor =	ncreader.get_dataset(step=step, name='y_vorticity')
z_vor =	ncreader.get_dataset(step=step, name='z_vorticity')

nz, ny, nx = x_vor.shape

vor = np.sqrt(x_vor ** 2 + y_vor ** 2 + z_vor ** 2)

print(t[step])


from mayavi import mlab
mlab.options.offscreen = True
mlab.figure(size=(2048, 1024))

#sf = mlab.pipeline.scalar_field(pres)
"""
obj = mlab.pipeline.volume(sf) #, vmin=0, vmax=0.8)
mlab.axes(xlabel=r'x', x_axis_visibility=True,
          ylabel=r'y', y_axis_visibility=True,
          zlabel=r'z', z_axis_visibility=True,
          ranges=[xmin, xmax, ymin, ymax, zmin, zmax])
mlab.colorbar(object=obj,
              title=r'pressure (m^2/s^2)',
              orientation='vertical')
mlab.savefig(filename='pressure.eps')
mlab.close()

mlab.figure(size=(1920, 1024))
obj = mlab.pipeline.iso_surface(sf, opacity=0.25)
mlab.axes(xlabel=r'x', x_axis_visibility=True,
          ylabel=r'y', y_axis_visibility=True,
          zlabel=r'z', z_axis_visibility=True,
          ranges=[xmin, xmax, ymin, ymax, zmin, zmax])
mlab.colorbar(object=obj,
              title=r'pressure (m^2/s^2)',
              orientation='vertical')

mlab.savefig(filename='iso_pressure.eps')
mlab.close()
"""

mlab.figure(size=(1024, 1024), bgcolor=(1,1,1), fgcolor=(0.,0.,0.))
sf = mlab.pipeline.scalar_field(vor)
#obj = mlab.pipeline.iso_surface(sf, opacity=0.4)
obj = mlab.pipeline.volume(sf)
ax = mlab.axes(xlabel='x', x_axis_visibility=True,
               ylabel='y', y_axis_visibility=True,
               zlabel='z', z_axis_visibility=True,
               ranges=[xmin, xmax, ymin, ymax, zmin, zmax])

ax.title_text_property.bold = False
ax.title_text_property.italic = True
ax.title_text_property.font_family = 'arial'
ax.label_text_property.bold = False
ax.label_text_property.italic = True
ax.label_text_property.font_family = 'arial'
ax.axes.font_factor = 1.0


colorbar = mlab.colorbar(object=obj,
              #title=r'vorticity magnitude (1/s)',
              label_fmt='%.2f',
              orientation='vertical')

# 1 July 2022
# https://github.com/enthought/mayavi/issues/643
colorbar.scalar_bar.unconstrained_font_size = True
colorbar.label_text_property.font_family = 'arial'
colorbar.label_text_property.bold = False
colorbar.label_text_property.font_size=24
colorbar.scalar_bar_representation.position = (0.9, 0.05)


mlab.view(azimuth=40.0, elevation=60, distance='auto', focalpoint='auto')
#print(mlab.view())

arr = mlab.screenshot()

fig = plt.figure(figsize=(16, 20), dpi=200)
grid = ImageGrid(fig, 111,
                 nrows_ncols=(3, 3),
                 axes_pad=(0.1, 0.24),
                 direction='column',
                 share_all=True,
                 cbar_location="bottom",
                 cbar_mode=None,
                 cbar_size="4%",
                 cbar_pad=0.8)

for i in range(9):
    ax = grid[i]
    ax.imshow(arr)
    ax.axis('off')
    if i < 3:
        add_timestamp(ax, t[step], xy=(0.01, 0.85), fmt="%.2f")

#plt.tight_layout()
plt.savefig('vorticity_magnitude.eps', bbox_inches='tight', dpi=200)

#mlab.savefig(filename='vorticity_magnitude.eps')
#mlab.close()

ncreader.close()
