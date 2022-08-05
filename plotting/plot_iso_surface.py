import numpy as np
import matplotlib.pyplot as plt
from nc_reader import nc_reader
from iso_surface import iso_surface
from mpl_toolkits.axes_grid1 import ImageGrid
from utils import *
import argparse
import os

parser = argparse.ArgumentParser(description='Create iso-surface plots.')
parser.add_argument('--filename',
                    type=str,
                    help='output file')
parser.add_argument('--step',
                    type=int)

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

parser.add_argument('--subfignum',
                    type=int,
                    help='subfigure number')

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
                    help='set vmin')

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


parser.add_argument('--color_points',
                    type=float,
                    nargs='+',
                    help='additional color points',
                    default=[])

parser.add_argument('--color_values',
                    type=float,
                    nargs='+',
                    help='rgb color values associated with the color points',
                    default=[])

parser.add_argument('--color_vmax',
                    type=float,
                    nargs=3,
                    help='rgb color at vmax',
                    default=None)

parser.add_argument('--color_vmin',
                    type=float,
                    nargs=3,
                    help='rgb color at vmin',
                    default=None)

parser.add_argument('--use_log_scale',
                    help='use log scale',
                    action='store_true')

parser.add_argument('--n_color_bar_ticks',
                    type=int,
                    default=10,
                    help='number of color bar ticks')


args = parser.parse_args()
fname = args.filename
step = args.step
enable_opacity = args.enable_opacity
invert_colormap = args.invert_colormap
save_path = args.save_path
overwrite = args.overwrite
fields = args.fields
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
n_color_bar_ticks = args.n_color_bar_ticks

if fields is None:
    print("No fields provided.")
    exit()


print()
print("\tFiles:           ", fname)
print("\tFields:          ", fields)
print("\tColor map:       ", colormap)
print("\tInvert color map:", invert_colormap)
print("\tEnable opacity:  ", enable_opacity)
print("\tOpacity vmax:    ", opacity_vmax)
print("\tOpacity vmin:    ", opacity_vmin)
print("\tColor vmax:      ", color_vmax)
print("\tColor vmin:      ", color_vmin)
print("\tColor points:    ", color_points)
print("\tColor values:    ", color_values)
print("\tColorbar vmin:   ", vmin)
print("\tColorbar ticks:  ", n_color_bar_ticks)
print("\tUse log scale:   ", use_log)
print("\tStep:            ", step)
print("\tSave path:       ", save_path)
print("\tOverwrite:       ", overwrite)
print("\tFignum:          ", fignum)
print()

ncreader = nc_reader()

for j, field in enumerate(fields):
    print('Plotting', field)

    ncreader.open(fname)
    t = ncreader.get_all('t')
    ncreader.close()

    iso = iso_surface(create_cmaps=True)
    iso.open(fname, add_time=False, width=1750, height=1600)
    iso.render(field_name=field, step=step,
               n_iso=100,
               vmin=vmin,
               colormap=colormap,
               enable_opacity=enable_opacity,
               opacity_vmax=opacity_vmax,
               opacity_vmin=opacity_vmin,
               opacity_points=opacity_points,
               opacity_values=opacity_values,
               invert_colormap=invert_colormap,
               color_vmax=color_vmax,
               color_vmin=color_vmin,
               color_points=color_points,
               color_values=color_values,
               use_log_scale=use_log,
               n_color_bar_ticks=n_color_bar_ticks,
               add_clabel=False)
    iso.export(file_path=save_path, file_name='temp_fig' + str(fignum) + '_' + str(subfignum) + '.png')
    iso.close()
    iso = None
