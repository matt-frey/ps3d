import numpy as np
import matplotlib.pyplot as plt
from nc_reader import nc_reader
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import ImageGrid
from utils import *
import argparse
import os

parser = argparse.ArgumentParser(description='Create iso-surface figure.')
parser.add_argument('--filename',
                    type=str,
                    help='output file')

parser.add_argument('--save_path',
                    type=str,
                    help='where to save the figures',
                    default=os.getcwd())

parser.add_argument('--steps',
                    type=int,
                    nargs=6,
                    help='list of steps',
                    default=[0, 1, 2, 3, 4, 5])

parser.add_argument('--n_iso',
                    type=int,
                    help='number of iso-surfaces',
                    default=20)

parser.add_argument('--overwrite',
                    help='overwrite figures',
                    action='store_true')

parser.add_argument('--fignum',
                    type=int,
                    help='figure number',
                    default=100)

args = parser.parse_args()
fname = args.filename
steps = args.steps
n_iso = args.n_iso
save_path = args.save_path
overwrite = args.overwrite
fignum = args.fignum

print()
print("\tFilename:               ", fname)
print("\tSteps:                  ", steps)
print("\tNumber of iso-surfaces: ", n_iso)
print("\tSave path:              ", save_path)
print("\tOverwrite:              ", overwrite)
print("\tFignum:                 ", fignum)
print()


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
    iso = iso_surface()

    figure = 'temp_figure.png'

    iso.open(fname)
    iso.render(step=step, n_iso=n_iso)
    iso.export(file_path=save_path, file_name=figure)
    iso.close()
    iso = None

    im = plt.imread(figure)

    grid[i].imshow(im)
    grid[i].axis('off')

    add_timestamp(grid[i], t[step], xy=(0.03, 1.06), fmt="%.2f")

save_figure(plt=plt, figpath=save_path, fignum=fignum, overwrite=overwrite)
plt.close()
