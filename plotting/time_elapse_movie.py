from iso_surface import iso_surface
import argparse
import os
from nc_reader import nc_reader

parser = argparse.ArgumentParser(description='Create evolution movie.')
parser.add_argument('--filenames',
                    type=str,
                    nargs='+'
                    help='NetCDF files')

parser.add_argument('--save_path',
                    type=str,
                    help='where to save PNGs and movie',
                    default=os.getcwd())

parser.add_argument('--field',
                    type=str,
                    help='which field to plot')

parser.add_argument('--starts',
                    type=int,
                    nargs='+',
                    help='start steps for each file')

parser.add_argument('--ends',
                    type=int,
                    nargs='+',
                    help='end steps for each file')

parser.add_argument('--movie_name',
                    type=str,
                    default='movie1.mp4',
                    help='name of movie file')

args = parser.parse_args()
fnames = args.filenames
field = args.field
start_step = args.start_step
end_step = args.end_step
save_path = args.save_path
movie_name = args.movie_name

print()
print("\tFilenames: ", fnames)
print("\tField:     ", field)
print("\tStart step:", start_step)
print("\tEnd step:  ", end_step)
print("\tSave path: ", save_path)
print("\tMovie name:", movie_name)
print()

basename = 'movie_fig_'

def render_step(fn, step, i):
    ncr = nc_reader()
    ncr.open(fn)
    t = ncr.get_all('t')[step]
    ncr.close()
    figure = basename + str(i).zfill(3) + '.png'
    iso = iso_surface(create_cmaps=False)
    iso.open(fname, width=1750, height=1600)
    iso.render(field_name=field, step=step,
               n_iso=100,
               vmin=0.0,
               colormap='rainbow4',
               enable_opacity=True,
               opacity_vmax=1.0,
               opacity_vmin=0.0,
               invert_colormap=True,
               n_color_bar_ticks=10,
               add_clabel=False)
    iso.export(file_path=save_path, file_name=figure)
    iso.close()
    del iso

    plt.figure(figsize=(8, 8), dpi=300)
    im = plt.imread(os.path.join(save_path, figure))
    plt.imshow(im)
    plt.axis('off')
    add_timestamp(plt, t, xy=(0.03, 1.06), fmt="%.2f")
    plt.savefig(os.path.join(save_path, figure), dpi=300, bbox_inches='tight')

def pngs2mp4:
    for fn in os.listdir(save_path):
        if 'movie_fig_' in fn:
            os.system('ffmpeg -i ' + os.path.join(save_path, basename + '%03d.png') +
                      ' -c:v libx264 -vf fps=25 ' + movie_name)


n = 0
for j, fname in enumerate(fnames):

    for step in range(start_step[j], end_step[j]+1):
        render_step(fname, step, n)

    n = n + 1

pngs2mp4()
