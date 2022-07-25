from tools.nc_reader import nc_reader
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
from utils import *
import argparse
import os

parser = argparse.ArgumentParser(description='Create power spectrum plot.')

parser.add_argument('--path',
                    type=str,
                    help='file path')

parser.add_argument('--ncfile',
                    type=str,
                    help='output file')

parser.add_argument('--efile',
                    type=str,
                    help='energy output file')

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
path = args.path
ncfile = args.ncfile
efile = args.efile
save_path = args.save_path
overwrite = args.overwrite
fignum = args.fignum

print()
print("\tPath:       ", path)
print("\tFilename:   ", ncfile)
print("\tEnergy file:", efile)
print("\tSave path:  ", save_path)
print("\tOverwrite:  ", overwrite)
print("\tFignum:     ", fignum)
print()


def find_steps:
    ts, kes, _ = np.loadtxt(os.path.join(fpath, efile), skiprows=1, unpack=True)
    kemax = kes.max()
    e1 = np.exp(1.0)
    e2 = np.exp(2.0)
    idx1 = find_nearest(kemax, kemax / e1)
    idx2 = find_nearest(kemax, kemax / e2)

    ncreader = nc_reader()
    ncreader.open(ncfile)
    t = ncreader.get_all('t')
    ncreader.close()

    step1 = find_nearest(t, ts[idx1])
    step2 = find_nearest(t, ts[idx2])

    return step1, step2

def create_spectra:
    step1, step2 = find_steps

    pathname = os.path.splitext(efile)[0]

    os.system('genspec --filename ' + efile + ' --step ' + str(step1))
    os.replace(os.path.join(path, pathname + '_spectrum.asc'),
               os.path.join(path, 'spectrum_exp1_decay.asc'))

    os.system('genspec --filename ' + efile + ' --step ' + str(step2))
    os.replace(os.path.join(path, pathname + '_spectrum.asc'),
               os.path.join(path, 'spectrum_exp2_decay.asc'))


def plot_spectrum(ax, ff):
    k, p = np.loadtxt(ff, skiprow=3 ,unpack=True)
    ax.loglog(k, p)

mpl.rcParams['font.size'] = 10

fig, axs = plt.subplots(2, 1, figsize=(8, 5), dpi=400, sharex=True, sharey=False)
grid = axs.flatten()

plot_spectrum(grid[0], os.path.join(path, 'spectrum_exp1_decay.asc'))
plot_spectrum(grid[1], os.path.join(path, 'spectrum_exp1_decay.asc'))

#for k in range(2):
    #grid[k].legend(loc='right', ncol=1, bbox_to_anchor=(1.28, 0.5))
#grid[1].set_xlabel(r'time, $t$')


plt.tight_layout()
save_figure(plt=plt, figpath=save_path, fignum=fignum, overwrite=overwrite)
plt.close()
