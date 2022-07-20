import numpy as np
import matplotlib.pyplot as plt
from nc_reader import nc_reader
import matplotlib as mpl
import colorcet as cc
from utils import *
import argparse
import os

parser = argparse.ArgumentParser(description='Volume rendering.')
parser.add_argument('--filename',
                    type=str,
                    help='output file')

parser.add_argument('--step',
                    type=int,
                    help='step number')

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
fname = args.filename
step = args.step
save_path = args.save_path
overwrite = args.overwrite
fignum = args.fignum

print()
print("\tFilename:  ", fname)
print("\tStep:      ", step)
print("\tSave path: ", save_path)
print("\tOverwrite: ", overwrite)
print("\tFignum:    ", fignum)
print()

mpl.rcParams['font.size'] = 16

ncreader = nc_reader()
ncreader.open(fname)

t = ncreader.get_all('t')

plt.figure(figsize=(9, 9))
image = make_volume_rendering(ncr=ncreader, step=step, field='vorticity_magnitude')
plt.imshow(image)
plt.axis('off')
add_timestamp(plt, t[step], xy=(0.05, 0.95), fmt="%.2f")

plt.tight_layout()

save_figure(plt=plt, figpath=save_path, fignum=fignum, overwrite=overwrite)
plt.close()

ncreader.close()
