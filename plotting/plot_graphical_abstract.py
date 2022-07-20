import numpy as np
import matplotlib.pyplot as plt
from nc_reader import nc_reader
import matplotlib as mpl
import colorcet as cc
from utils import *
import argparse
import os

parser = argparse.ArgumentParser(description='Create graphical abstract.')
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

args = parser.parse_args()
fname = args.filename
step = args.step
save_path = args.save_path

print()
print("\tFilename:  ", fname)
print("\tStep:      ", step)
print("\tSave path: ", save_path)
print()

mpl.rcParams['font.size'] = 16

ncreader = nc_reader()
ncreader.open(fname)

#t = ncreader.get_all('t')

# ratio 1.2:1
# 1 cm = 37.7952755906 pixel
px = 37.7952755906
height_pixels = 5 * px
image = make_volume_rendering(ncr=ncreader, step=step, field='vorticity_magnitude',
                              show_axes=False, show_colorbar=False, fmt='jpeg',
                              surface_count=150, width=1.2*height_pixels, height=height_pixels,
                              scale=10.0, margin={'r': 1, 'b': 1, 'l': 1, 't': 1})

os.rename('temp_figure.jpeg', 'graphical_abstract.jpeg')
plt.close()

ncreader.close()
