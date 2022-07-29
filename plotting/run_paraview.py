from iso_surface import iso_surface
import argparse
import os

parser = argparse.ArgumentParser(description='Create iso-surface figures.')
parser.add_argument('--filename',
                    type=str,
                    help='output file')

parser.add_argument('--save_path',
                    type=str,
                    help='where to save the figures',
                    default=os.getcwd())

parser.add_argument('--steps',
                    type=int,
                    nargs='+',
                    help='list of steps',
                    default=[0, 1])

parser.add_argument('--n_iso',
                    type=int,
                    help='number of iso-surfaces',
                    default=20)

parser.add_argument('--overwrite',
                    help='overwrite figures',
                    action='store_true')

parser.add_argument('--fignums',
                    type=int,
                    nargs='+',
                    help='figure numbers',
                    default=[100, 101])

args = parser.parse_args()
fname = args.filename
steps = args.steps
n_iso = args.n_iso
save_path = args.save_path
overwrite = args.overwrite
fignums = args.fignums

print()
print("\tFilename:               ", fname)
print("\tSteps:                  ", steps)
print("\tNumber of iso-surfaces: ", n_iso)
print("\tSave path:              ", save_path)
print("\tOverwrite:              ", overwrite)
print("\tFignums:                ", fignums)
print()

if not len(fignums) == len(steps):
    print("Error: The number of figures and the number of steps not identical.")


for i, step in enumerate(steps):
    iso = iso_surface()

    figure = 'fig' + str(fignums[i]) + '.png'
    full_path = os.path.join(save_path, figure)

    if os.path.exists(full_path) and overwrite:
        print("File '" + full_path + "' already exists. Exiting.")
        exit()

    iso.open(fname)
    iso.render(step=step, n_iso=n_iso)
    iso.export(file_path=save_path, file_name='fig' + str(fignums[i]) + '.png')
    iso.close()
    iso = None

#iso.save_camera_orbiting_animation(step=60, n_frames=360, file_path='./',
                                   #file_name="test.mp4", keep_frames=False)


#iso.save_animation(beg=0, end=100,
                   #fps=25,
                   #file_path='./',
                   #file_name="animation.mp4",
                   #keep_frames=True)

