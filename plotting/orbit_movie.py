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
step = args.step
save_path = args.save_path

print()
print("\tFilename:               ", fname)
print("\tField:                  ", field)
print("\tStep:                   ", step)
print("\tSave path:              ", save_path)
print()

iso = iso_surface(create_cmaps=True)
iso.open(fname, width=1750, height=1600)

iso.save_camera_orbiting_animation(field_name=field,
                                   step=step,
                                   n_frames=360,
                                   file_path='./movies',
                                   file_name="movie1.mp4",
                                   keep_frames=False,
                                   n_iso=100,
                                   vmin=0.0,
                                   colormap='rainbow4',
                                   enable_opacity=True,
                                   opacity_vmax=1.0,
                                   opacity_vmin=0.0,
                                   invert_colormap=True)

