from iso_surface import iso_surface
import argparse
import os

parser = argparse.ArgumentParser(description='Create orbiting movie.')
parser.add_argument('--filename',
                    type=str,
                    help='output file')

parser.add_argument('--save_path',
                    type=str,
                    help='where to save the figures',
                    default=os.getcwd())

parser.add_argument('--field',
                    type=str,
                    help='which field to plot')

parser.add_argument('--step',
                    type=int,
                    help='step number',
                    default=0)

args = parser.parse_args()
fname = args.filename
field = args.field
steps = args.steps
save_path = args.save_path

print()
print("\tFilename:               ", fname)
print("\tField:                  ", field)
print("\tStep:                   ", step)
print("\tSave path:              ", save_path)
print()

iso.save_camera_orbiting_animation(step=step, n_frames=360,
                                   file_path='./movies',
                                   file_name="movie1.mp4",
                                   keep_frames=False)


#iso.save_animation(beg=0, end=100,
                   #fps=25,
                   #file_path='./',
                   #file_name="animation.mp4",
                   #keep_frames=True)

