import numpy as np
import matplotlib.pyplot as plt
from nc_reader import nc_reader
import matplotlib as mpl
import colorcet as cc
from utils import *
import argparse
import os
from iso_surface import *

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

iso = iso_surface(create_cmaps=True)

# ratio must be 1.2 : 1
iso.open(fname, add_time=False, width=1920, height=1600,
         add_axes=False)
iso.render(field_name='vorticity_magnitude', step=step,
           n_iso=100,
           vmin=0.0,
           colormap='rainbow4',
           enable_opacity=True,
           opacity_vmax=1.0,
           opacity_vmin=0.0,
           invert_colormap=True,
           cam_focal_point=[0, 0, 0],
           add_color_bar=False)
iso.export(file_path=save_path, file_name='graphical_abstract.jpeg', quality=100)
iso.close()
iso = None
