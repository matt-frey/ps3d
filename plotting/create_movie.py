import numpy as np
import matplotlib.pyplot as plt
from nc_reader import nc_reader
import matplotlib as mpl
import colorcet as cc
from utils import *
import argparse
import os

parser = argparse.ArgumentParser(description='Create movie.')
parser.add_argument('--path',
                    type=str,
                    help='path to files')

parser.add_argument('--save_path',
                    type=str,
                    help='where to save the figures',
                    default=os.getcwd())

parser.add_argument('--overwrite',
                    help='overwrite figures',
                    action='store_true')

args = parser.parse_args()
path = args.path
save_path = args.save_path
overwrite = args.overwrite

print()
print("\tPath:      ", path)
print("\tSave path: ", save_path)
print("\tOverwrite: ", overwrite)
print()

i = 0

def create_images(fn, lo, hi):
    ncreader = nc_reader()
    ncreader.open(fn)
    t = ncreader.get_all('t')
    for step in range(lo, hi)
        print("Creating plot", i, "for time", t[step], "...", end='')
        plt.figure(figsize=(9, 9))
        image = make_volume_rendering(ncr=ncreader, step=step, field='vorticity_magnitude', fmt='png')
        plt.imshow(image)
        plt.axis('off')
        add_timestamp(plt, t[step], xy=(0.05, 0.95), fmt="%.2f")
        plt.tight_layout()
        fname = 'fig' + str(i).zfill(5) + '.jpeg'
        plt.savefig(fname=fname, format='jpeg')
        plt.close()
        os.replace(fname, os.path.join('movie_temp_dir', fname))
        i = i + 1
        print("done.")
    ncreader.close()

mpl.rcParams['font.size'] = 16

fname1 'beltrami_256_fields.nc'
fname2 = 'beltrami_256_restart_fields.nc'

if not os.path.exists('movie_temp_dir'):
    os.mkdir('movie_temp_dir')



create_images(fname1, 0, 5)
create_images(fname2, 0, 41)
create_images(fname1, 7, 11)


#sudo apt install ffmpeg
#ffmpeg -r 1/5 -i fig%05d.jpeg -c:v libx264 -vf fps=25 -pix_fmt yuv420p out.mp4

ncreader.close()
