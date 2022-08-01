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
                    nargs=3,
                    default=[0, 1, 2],
                    help='add 3 steps to plot')

parser.add_argument('--file_numbers',
                    type=int,
                    nargs=3,
                    default=[0, 0, 0],
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

args = parser.parse_args()
fnames = args.filenames
steps = np.asarray(args.steps)
file_numbers = np.asarray(args.file_numbers)
save_path = args.save_path
overwrite = args.overwrite
fignum = args.fignum

print()
print("\tFiles:           ", fnames)
print("\tSteps:           ", steps)
print("\tFile numbers:    ", file_numbers)
print("\tSave path:       ", save_path)
print("\tOverwrite:       ", overwrite)
print("\tFignum:          ", fignum)
print()

ncreader = nc_reader()

fig = plt.figure(figsize=(7.5, 7), dpi=500)
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

# create pressure images:
for i in range(3):
    step = steps[i]
    
    iso = iso_surface()
    iso.open(fnames[file_numbers[i]], add_time=False, width=1750, height=1600)
    iso.render(field_name='pressure', step=step,
               n_iso=40,
               colormap='Cool to Warm',
               enable_opacity=True,
               opacity_vmax=1.0,
               opacity_vmin=1.0,
               opacity_points=[0.0],
               opacity_values=[0.0])
    iso.export(file_path=save_path, file_name='temp_fig' + str(fignum) + '_' + str(i) + '.png')
    iso.close()
    iso = None

#create helicity figures
for i in range(3):
    step = steps[i]

    iso = iso_surface()
    iso.open(fnames[file_numbers[i]], add_time=False, width=1750, height=1600)
    iso.render(field_name='helicity', step=step,
               n_iso=40,
               colormap='Cool to warm',
               enable_opacity=True,
               opacity_vmax=1.0,
               opacity_vmin=1.0,
               opacity_points=[0.0],
               opacity_values=[0.0])
    iso.export(file_path=save_path, file_name='temp_fig' + str(fignum) + '_' + str(i+3) + '.png')
    iso.close()
    iso = None


for i in range(3):
    ncreader.open(fnames[file_numbers[i]])
    t = ncreader.get_all('t')
    nreader.close()
    step = steps[i]
    
    im = plt.imread(os.path.join(save_path, 'temp_fig' + str(fignum) + '_' + str(i) + '.png'))
    ax = grid[i]
    ax.imshow(im)
    ax.axis('off')
    add_timestamp(ax, t[step], xy=(0.2, 0.95), fmt="%.2f", fontsize=8)

for i in range(3):
    ncreader.open(fnames[file_numbers[i]])
    t = ncreader.get_all('t')
    nreader.close()
    step = steps[i]

    im = plt.imread(os.path.join(save_path, 'temp_fig' + str(fignum) + '_' + str(i+3) + '.png'))
    ax = grid[i+3]
    ax.imshow(im)
    ax.axis('off')
    add_timestamp(ax, t[step], xy=(0.2, 0.95), fmt="%.2f", fontsize=8)
   
save_figure(plt=plt, figpath=save_path, fignum=fignum, overwrite=overwrite)
plt.close()

for i in range(6):
    os.remove(os.path.join(save_path, 'temp_fig' + str(fignum) + '_' + str(i) + '.png'))
