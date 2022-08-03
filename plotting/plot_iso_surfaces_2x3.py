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
                    nargs=6,
                    default=[0, 1, 2, 3, 4, 5],
                    help='add 3 steps to plot')

parser.add_argument('--file_numbers',
                    type=int,
                    nargs=6,
                    default=[0, 0, 0, 0, 0, 0],
                    help='the file by index to read data from')

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

parser.add_argument('--colormap',
                    default='Cool to Warm',
                    help='color map')

parser.add_argument('--enable_opacity',
                    help='activate opacity',
                    action='store_true')

parser.add_argument('--vmin',
                    type=str,
                    help='set vmin')

parser.add_argument('--invert_colormap',
                    help='invert the colormap',
                    action='store_true')

parser.add_argument('--opacity_points',
                    type=str,
                    help='additional opacity points w.r.t. vmax',
                    default=[])

parser.add_argument('--opacity_values',
                    type=str,
                    help='opacity values associated with the opacity points',
                    default=[])

parser.add_argument('--opacity_vmax',
                    type=str,
                    help='opacity at vmax')

parser.add_argument('--opacity_vmin',
                    type=str,
                    help='opacity at vmin')

parser.add_argument('--field',
                    help='field to plot')

parser.add_argument('--color_points',
                    type=str,
                    nargs='+',
                    help='additional color points',
                    default=[])

parser.add_argument('--color_values',
                    type=str,
                    nargs='+',
                    help='rgb color values associated with the color points',
                    default=[])

parser.add_argument('--color_vmax',
                    type=str,
                    help='rgb color at vmax',
                    default=None)

parser.add_argument('--color_vmin',
                    type=str,
                    help='rgb color at vmin',
                    default=None)

parser.add_argument('--use_log_scale',
                    nargs=6,
                    help='use log scale',
                    action='store_true')

args = parser.parse_args()
fnames = args.filenames
steps = np.asarray(args.steps)
enable_opacity = ""
if args.enable_opacity:
    enable_opacity = ' --enable_opacity'

invert_colormap = ''
if args.invert_colormap:
    invert_colormap = ' --invert_colormap'

save_path = args.save_path
overwrite = args.overwrite
field = args.field
fignum = args.fignum
subfignum = args.subfignum
colormap = args.colormap
use_log = args.use_log_scale
vmin = args.vmin
opacity_points = args.opacity_points
opacity_values = args.opacity_values
opacity_vmax = args.opacity_vmax
opacity_vmin = args.opacity_vmin
color_points = args.color_points
color_values = args.color_values
color_vmax = args.color_vmax
color_vmin = args.color_vmin
file_numbers = np.asarray(args.file_numbers)
save_path = args.save_path
overwrite = args.overwrite
fignum = args.fignum

print()
print("\tFiles:           ", fnames)
print("\tField:           ", field)
print("\tColor map:       ", colormap)
print("\tInvert color map:", args.invert_colormap)
print("\tEnable opacity:  ", args.senable_opacity)
print("\tOpacity vmax:    ", opacity_vmax)
print("\tOpacity vmin:    ", opacity_vmin)
print("\tColor vmax:      ", color_vmax)
print("\tColor vmin:      ", color_vmin)
print("\tColor points:    ", color_points)
print("\tColor values:    ", color_values)
print("\tColorbar vmin:   ", vmin)
print("\tUse log scale:   ", use_log)
print("\tSteps:           ", steps)
print("\tFile numbers:    ", file_numbers)
print("\tSave path:       ", save_path)
print("\tOverwrite:       ", overwrite)
print("\tFignum:          ", fignum)
print()



for i in range(6):
    do_log = ''
    if use_log[i]:
        do_log = ' --use_log_scale'
    os.system("python " + os.path.join(spath, "plot_iso_surface.py") + \
              " --filename " + os.path.join(fpath, fnames[file_numbers[i]]) + \
              " --step " + str(steps[i]) + \
              " --subfignum " + str(i) + \
              " --field " + field + \
              " --colormap " + colormap + \
              + invert_colormap + \
              + enable_opacity + \
              " --vmin " + vmin + \
              " --opacity_points " + opacity_points + \
              " --opacity_values " + opacity_values + \
              " --opacity_vmax " + opacity_vmax + \
              " --opacity_vmin " + opacity_vmax + \
              " --color_points " + color_points + \
              " --color_values " + color_values + \
              " --color_vmax " + color_vmax + \
              " --color_vmin " + color_vmin + \
              + do_log + \
              " --fignum " + str(fignum) + \
              " --overwrite" + \
              " --save_path " + save_path)

ncreader = nc_reader()

fig = plt.figure(figsize=(7.5, 4.5), dpi=500)
grid = ImageGrid(fig, 111,
                 nrows_ncols=(2, 3),
                 aspect=True,
                 axes_pad=(0.01, 0.01),
                 direction='row',
                 share_all=True,
                 cbar_location="right",
                 cbar_mode='none',
                 cbar_size="4%",
                 cbar_pad=0.05)

for i in range(6):
    ncreader.open(fnames[file_numbers[i]])
    t = ncreader.get_all('t')
    ncreader.close()
    step = steps[i]

    im = plt.imread(os.path.join(save_path, 'temp_fig' + str(fignum) + '_' + str(i) + '.png'))
    ax = grid[i]
    ax.imshow(im)
    ax.axis('off')
    add_timestamp(ax, t[step], xy=(0.2, 0.95), fmt="%.2f", fontsize=8)

save_figure(plt=plt, figpath=save_path, fignum=fignum, overwrite=overwrite)
plt.close()

#for i in range(6):
#    os.remove(os.path.join(save_path, 'temp_fig' + str(fignum) + '_' + str(i) + '.png'))
