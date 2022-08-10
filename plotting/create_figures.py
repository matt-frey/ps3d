#python create_figures.py --file_path paper_runs/ --save_path ./figures/ --script_path $HOME/ps3d/plotting
import argparse
import numpy as np
import os

parser = argparse.ArgumentParser(
    description='Create all plots.')

parser.add_argument('--file_path',
                    type=str,
                    help='path to files')

parser.add_argument('--script_path',
                    type=str,
                    help='path to python scripts')

parser.add_argument('--save_path',
                    type=str,
                    help='where to save the figures',
                    default=os.getcwd())

parser.add_argument('--overwrite',
                    help='overwrite figures',
                    action='store_true')

parser.add_argument('--figures',
                    nargs='*',
                    type=int,
                    help='which figures to plot',
                    default=np.arange(2, 16, dtype=int))

parser.add_argument('--movies',
                    help='create movies',
                    action='store_true')

parser.add_argument('--graphical_abstract',
                    help='create graphical abstract',
                    action='store_true')

args = parser.parse_args()
fpath = args.file_path
spath = args.script_path
save_path = args.save_path
overwrite = args.overwrite
figures = args.figures
movies = args.movies
ga = args.graphical_abstract

# Figure 2:
if 2 in figures:
    os.system("python " + os.path.join(spath, "plot_ke_en_varying_prediss.py") + \
              " --filepath " + fpath + \
              " --fignum 2" + \
              " --overwrite" + \
              " --save_path " + save_path)

# Figure 3:
if 3 in figures:
    os.system("python " + os.path.join(spath, "plot_grid_resolution.py") + \
              " --filepath " + fpath + \
              " --fignum 3" + \
              " --overwrite" + \
              " --save_path " + save_path)

# Figure 4:
if 4 in figures:
    os.system("python " + os.path.join(spath, "plot_slice_difference.py") + \
              " --filename " + os.path.join(fpath, "beltrami_256_fields.nc") + \
              " --step 5" + \
              " --plane 'xz'" + \
              " --loc 128" + \
              " --zlabel '$y = 0$'" + \
              " --fignum 4" +  \
              " --overwrite" + \
              " --save_path " + save_path)

# Figure 5:
if 5 in figures:
    os.system("python " + os.path.join(spath, "plot_slice_difference.py") + \
              " --filename " + os.path.join(fpath, "beltrami_256_fields.nc") + \
              " --step 5" + \
              " --plane 'xy'" + \
              " --loc 128" + \
              " --zlabel '$z = 0$'" + \
            " --fignum 5" + \
              " --overwrite" + \
              " --save_path " + save_path)

# Figure 6:
if 6 in figures:
    #>>> t, ke, en = np.loadtxt('paper_runs/beltrami_256_restart_ecomp.asc', unpack=True)
    #>>> idx = np.argmax(en)
    #>>> t[idx]
    #62.916624
    # t = 63 is step 26
    for i, step in enumerate([0, 10, 16, 20, 26, 40]):
        os.system("python " + os.path.join(spath, "plot_iso_surface.py") + \
                  " --filename " + os.path.join(fpath, "beltrami_256_restart_fields.nc") + \
                  " --step " + str(step) + \
                  " --subfignum " + str(i) + \
                  " --fields vorticity_magnitude" + \
                  " --colormap rainbow4" + \
                  " --invert_colormap" + \
                  " --enable_opacity" + \
                  " --vmin 0.0" + \
                  " --opacity_vmax 1.0" + \
                  " --opacity_vmin 0.0" + \
                  " --fignum 6" + \
                  " --overwrite" + \
                  " --save_path " + save_path)


    os.system("python " + os.path.join(spath, "plot_iso_surfaces_2x3.py") + \
              " --filenames " + os.path.join(fpath, "beltrami_256_restart_fields.nc ") + \
              " --steps 0 10 16 20 26 40" + \
              " --file_numbers 0 0 0 0 0 0" + \
              " --fignum 6" + \
              " --overwrite" + \
              " --save_path " + save_path)

# Figure 7:
# create xy-plane cross-section of vorticty magnitude
# at iy = 128 (i.e. y = 0):
if 7 in figures:
    os.system("python " + os.path.join(spath, "plot_slices_2x3.py") + \
              " --filenames " + os.path.join(fpath, "beltrami_256_restart_fields.nc") + \
              " --plane 'xz'" + \
              " --loc 128" + \
              " --zlabel '$y = 0$'" + \
              " --steps 0 10 16 20 26 40" + \
              " --file_numbers 0 0 0 0 0 0 " + \
              " --field vorticity_magnitude" + \
              " --colormaps rainbow4_r rainbow4_r rainbow4_r rainbow4_r rainbow4_r rainbow4_r" + \
              " --norms none none none log log log" + \
              " --fignum 7" + \
              " --overwrite" + \
              " --save_path " + save_path)

# Figure 8:
if 8 in figures:
    os.system("python " + os.path.join(spath, "plot_slices_2x3.py") + \
              " --filenames " + os.path.join(fpath, "beltrami_256_restart_fields.nc") + \
              " --plane 'xy'" + \
              " --loc 128" + \
              " --zlabel '$z = 0$'" + \
              " --steps 0 10 16 20 26 40" + \
              " --file_numbers 0 0 0 0 0 0" + \
              " --field vorticity_magnitude" + \
              " --colormaps rainbow4_r rainbow4_r rainbow4_r rainbow4_r rainbow4_r rainbow4_r" + \
              " --norms none none none log log log" + \
              " --fignum 8" + \
              " --overwrite" + \
              " --save_path " + save_path)

# Figure 9:
if 9 in figures:
    os.system("python " + os.path.join(spath, "plot_slices_2x3.py") + \
              " --filenames " + os.path.join(fpath, "beltrami_256_restart_fields.nc") + \
              " --plane 'xy'" + \
              " --loc 256" + \
              " --zlabel '$z = \pi/2$'" + \
              " --steps 0 10 16 20 26 40" + \
              " --file_numbers 0 0 0 0 0 0" + \
              " --field vorticity_magnitude" + \
              " --colormaps rainbow4_r rainbow4_r rainbow4_r rainbow4_r rainbow4_r rainbow4_r" + \
              " --norms none none none log log log" + \
              " --fignum 9" + \
              " --overwrite" + \
              " --save_path " + save_path)


# Figure 10:
if 10 in figures:
    for i, step in enumerate([0, 10, 16, 20, 26, 40]):
        os.system("python " + os.path.join(spath, "plot_iso_surface.py") + \
                  " --filename " + os.path.join(fpath, "beltrami_256_restart_fields.nc") + \
                  " --step " + str(step) + \
                  " --subfignum " + str(i) + \
                  " --fields helicity" + \
                  " --colormap rainbow4" + \
                  " --invert_colormap" + \
                  " --enable_opacity" + \
                  " --opacity_vmax 1.0" + \
                  " --opacity_vmin 1.0" + \
	          " --opacity_points 0.0" +	\
	          " --opacity_values 0.0" +	\
                  " --fignum 10" + \
                  " --overwrite" + \
                  " --save_path " + save_path)

    os.system("python " + os.path.join(spath, "plot_iso_surfaces_2x3.py") + \
              " --filenames " + os.path.join(fpath, "beltrami_256_restart_fields.nc") + \
              " --steps 0 10 16 20 26 40" + \
              " --file_numbers 0 0 0 0 0 0" + \
              " --fignum 10" + \
              " --overwrite" + \
              " --save_path " + save_path)

# Figure 11:
if 11 in figures:
    os.system("python " + os.path.join(spath, "plot_slices_2x3.py") + \
              " --filenames " + os.path.join(fpath, "beltrami_256_restart_fields.nc") + \
              " --plane 'xy'" + \
              " --loc 256" + \
              " --zlabel '$z = \pi/2$'" + \
              " --steps 0 10 16 20 26 40" + \
              " --file_numbers 0 0 0 0 0 0" + \
              " --field helicity" + \
              " --colormaps rainbow4_r rainbow4_r rainbow4_r rainbow4_r rainbow4_r rainbow4_r" + \
              " --norms centered centered centered centered centered centered" + \
              " --fignum 11" + \
              " --overwrite" + \
              " --save_path " + save_path)

# Figure 12:
if 12 in figures:
    os.system("python " + os.path.join(spath, "plot_vor_vel_profiles.py") + \
              " --filenames " + os.path.join(fpath, "beltrami_256_restart_fields.nc ") + \
              os.path.join(fpath, "beltrami_256_fields.nc") + \
              " --steps 0 10 16 20 26 40 10" + \
              " --file_numbers 0 0 0 0 0 0 1" + \
              " --fignum 12" + \
              " --overwrite" + \
              " --save_path " + save_path)

# Figure 13:
# create velocity, vorticity and helicity evolution:
if 13 in figures:
    os.system("python " + os.path.join(spath, "plot_vor_vel_he_evolution.py") + \
              " --filename " + os.path.join(fpath, "beltrami_256_fields.nc") + \
              " --restartfile " + os.path.join(fpath, "beltrami_256_restart_fields.nc") + \
              " --fignum 110 13" + \
              " --overwrite" + \
              " --save_path " + save_path)

# Figure 14:
# create power spectrum plot
if 14 in figures:
    os.system("python " + os.path.join(spath, "plot_power_spectra.py") + \
              " --path " + fpath + \
              " --ncfile beltrami_256_fields.nc beltrami_256_restart_fields.nc" + \
              " --efile  beltrami_256_ecomp.asc" + \
              " --fignum 14" + \
              " --overwrite" + \
              " --save_path " + save_path)

if 16 in figures:
    os.system("python " + os.path.join(spath, "plot_vor2vel_conv.py") + \
              " --path " + os.path.join(fpath, "..", "tests") + \
              " --fignum 16" + \
              " --overwrite" + \
              " --save_path " + save_path)

if 17 in figures:
    os.system("python " + os.path.join(spath, "plot_filter_hyper.py") + \
              " --path " + os.path.join(fpath, "..", "tests") + \
              " --labels '$nz = 128$' '$nz = 256$'" + \
              " --fignum 170" + \
              " --overwrite" + \
              " --save_path " + save_path)


if 150 in figures:
    os.system("python " + os.path.join(spath, "plot_contours.py") + \
              " --filename " + os.path.join(fpath, "beltrami_256_restart_fields.nc") + \
              " --step 26"
              " --fignum 150" + \
              " --overwrite" + \
              " --save_path " + save_path)

if 130 in figures:
    os.system("python " + os.path.join(spath, "plot_volume_fraction.py") + \
              " --filename " + os.path.join(fpath, "beltrami_256_fields.nc") + \
              " --restartfile " + os.path.join(fpath, "beltrami_256_restart_fields.nc") + \
              " --fignum 130" + \
              " --overwrite" + \
              " --save_path " + save_path)


if 100 in figures:
    os.system("python " + os.path.join(spath, "plot_diff_velocity_rms_evolution.py") + \
              " --filename " + os.path.join(fpath, "beltrami_256_restart_fields.nc") + \
              " --fignum 100" + \
              " --overwrite" + \
              " --save_path " + save_path)

if movies:
    os.system("python " + os.path.join(spath, "orbit_movie.py") + \
              " --filename " + os.path.join(fpath, "beltrami_256_restart_fields.nc") + \
              " --step 26" + \
              " --field vorticity_magnitude" + \
              " --save_path " + save_path)

if ga:
    os.system("python " + os.path.join(spath, "plot_graphical_abstract.py") + \
              " --filename " + os.path.join(fpath, "beltrami_256_restart_fields.nc") + \
              " --step 20" + \
              " --save_path " + save_path)
