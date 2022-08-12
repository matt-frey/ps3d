import argparse
import os
from nc_reader import nc_reader
from utils import *
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description='Create evolution movie.')
parser.add_argument('--filenames',
                    type=str,
                    nargs='+',
                    help='NetCDF files')

parser.add_argument('--save_path',
                    type=str,
                    help='where to save PNGs and movie',
                    default=os.getcwd())

parser.add_argument('--field',
                    type=str,
                    help='which field to plot')

parser.add_argument('--start_step',
                    type=int,
                    nargs='+',
                    help='start steps for each file')

parser.add_argument('--end_step',
                    type=int,
                    nargs='+',
                    help='end steps for each file')

parser.add_argument('--movie_name',
                    type=str,
                    default='movie1.mp4',
                    help='name of movie file')
parser.add_argument('--script_path',
                    type=str,
                    help='Path to plot_iso_surface.py')

args = parser.parse_args()
fnames = args.filenames
field = args.field
start_step = args.start_step
end_step = args.end_step
save_path = args.save_path
movie_name = args.movie_name
spath = args.script_path

print()
print("\tFilenames: ", fnames)
print("\tField:     ", field)
print("\tStart step:", start_step)
print("\tEnd step:  ", end_step)
print("\tSave path: ", save_path)
print("\tMovie name:", movie_name)
print()

fignum = 1000
basename = 'temp_fig' + str(fignum)

def render_step(fn, step, i):
     os.system("python " + os.path.join(spath, "plot_iso_surface.py") + \
              " --filename " + fn + \
              " --step " + str(step) + \
              " --subfignum " + str(i) + \
              " --fields " + field + \
              " --colormap rainbow4" + \
              " --invert_colormap" + \
              " --enable_opacity" + \
              " --vmin 0.0" + \
              " --opacity_vmax 1.0" + \
              " --opacity_vmin 0.0" + \
              " --fignum " + str(fignum) + \
              " --overwrite" + \
              " --add_clabel True" + \
              " --add_time True" + \
              " --time_format 4.2f" + \
              " --save_path " + save_path)
     os.rename(os.path.join(save_path, basename + '_' + str(i) + '.png'),
               os.path.join(save_path, basename + '_' + str(i).zfill(3) + '.png'))

n = 0
for j, fname in enumerate(fnames):
    for step in range(start_step[j], end_step[j]+1):
        render_step(fname, step, n)
        n = n + 1

os.system('ffmpeg -r 1.5 -i ' + os.path.join(save_path, basename + '_%03d.png') + \
          ' -c:v libx264 -vf fps=25 ' + movie_name)
