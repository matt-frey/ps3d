import numpy as np
import matplotlib.pyplot as plt
from nc_reader import nc_reader
from iso_surface import iso_surface
from mpl_toolkits.axes_grid1 import ImageGrid
from utils import *
import argparse
import os

parser = argparse.ArgumentParser(description='Create iso-surface plots.')
parser.add_argument('--filenames',
                    type=str,
                    nargs='+',
                    help='output file')
parser.add_argument('--steps',
                    type=int,
                    nargs=9,
                    default=[0, 1, 2, 3, 4, 5, 6, 7, 8],
                    help='add 9 steps to plot')

parser.add_argument('--file_numbers',
                    type=int,
                    nargs=9,
                    default=[0, 0, 0, 0, 0, 0, 0, 0, 0],
                    help='the file by index to read data from')

parser.add_argument('--save_path',
                    type=str,
                    help='where to save the figures',
                    default=os.getcwd())

parser.add_argument('--overwrite',
                    help='overwrite figures',
                    action='store_true')

parser.add_argument('--fignums',
                    type=int,
                    nargs='+',
                    help='figure numbers')

parser.add_argument('--fields',
                    nargs='+',
                    help='fields to plot')

parser.add_argument('--colormap',
                    default='Cool to Warm',
                    help='color map')

parser.add_argument('--enable_opacity',
                    help='activate opacity',
                    action='store_true')

parser.add_argument('--vmin',
                    type=float,
                    help='set vmin',
                    default=None)

parser.add_argument('--invert_colormap',
                    help='invert the colormap',
                    action='store_true')

parser.add_argument('--opacity_points',
                    type=float,
                    nargs='+',
                    help='additional opacity points w.r.t. vmax',
                    default=[])

parser.add_argument('--opacity_values',
                    type=float,
                    nargs='+',
                    help='opacity values associated with the opacity points',
                    default=[])

parser.add_argument('--opacity_vmax',
                    type=float,
                    help='opacity at vmax')

parser.add_argument('--opacity_vmin',
                    type=float,
                    help='opacity at vmin')

args = parser.parse_args()
fnames = args.filenames
steps = np.asarray(args.steps)
file_numbers = np.asarray(args.file_numbers)
enable_opacity = args.enable_opacity
invert_colormap = args.invert_colormap
save_path = args.save_path
overwrite = args.overwrite
fields = args.fields
fignums = args.fignums
colormap = args.colormap
vmin = args.vmin
opacity_points = args.opacity_points
opacity_values = args.opacity_values
opacity_vmax = args.opacity_vmax
opacity_vmin = args.opacity_vmin

if fields is None:
    print("No fields provided.")
    exit()

if fignums is None:
    fignums = np.arange(1, len(fields) + 1, dtype=int)

if not len(fields) == len(fignums):
    print("Number of fields and number of figures do not match.")
    exit()

print()
print("\tFiles:           ", fnames)
print("\tFields:          ", fields)
print("\tColor map:       ", colormap)
print("\tInvert color map:", invert_colormap)
print("\tEnable opacity:  ", enable_opacity)
print("\tOpacity vmax:    ", opacity_vmax)
print("\tOpacity vmin:    ", opacity_vmin)
print("\tColorbar vmin:   ", vmin)
print("\tSteps:           ", steps)
print("\tFile numbers:    ", file_numbers)
print("\tSave path:       ", save_path)
print("\tOverwrite:       ", overwrite)
print("\tFignums:         ", fignums)
print()

ncreader = nc_reader()

for j, field in enumerate(fields):
    print('Plotting', field)
    fig = plt.figure(figsize=(7.5, 7), dpi=500)
    grid = ImageGrid(fig, 111,
                    nrows_ncols=(3, 3),
                    aspect=True,
                    axes_pad=(0.01, 0.01),
                    direction='row',
                    share_all=True,
                    cbar_location="right",
                    cbar_mode='none',
                    cbar_size="4%",
                    cbar_pad=0.05)

    for i, step in enumerate(steps):
        ncreader.open(fnames[file_numbers[i]])
        t = ncreader.get_all('t')
        ncreader.close()

        iso = iso_surface()
        iso.open(fnames[file_numbers[i]], add_time=False, width=1750, height=1600)
        iso.render(field_name=field, step=step,
                   n_iso=40,
                   vmin=vmin,
                   colormap=colormap,
                   enable_opacity=enable_opacity,
                   opacity_vmax=opacity_vmax,
                   opacity_vmin=opacity_vmin,
                   opacity_points=opacity_points,
                   opacity_values=opacity_values,
                   invert_colormap=invert_colormap)
        iso.export(file_path=save_path, file_name='temp_fig' + str(fignums[j]) + '_' + str(i) + '.png')
        iso.close()
        iso = None

    for i, step in enumerate(steps):
        im = plt.imread(os.path.join(save_path, 'temp_fig' + str(fignums[j]) + '_' + str(i) + '.png'))
        ax = grid[i]
        ax.imshow(im)
        ax.axis('off')
        add_timestamp(ax, t[step], xy=(0.2, 0.95), fmt="%.2f", fontsize=8)

    save_figure(plt=plt, figpath=save_path, fignum=fignums[j], overwrite=overwrite)
    plt.close()

    for i in range(9):
        os.remove(os.path.join(save_path, 'temp_fig' + str(fignums[j]) + '_' + str(i) + '.png'))
