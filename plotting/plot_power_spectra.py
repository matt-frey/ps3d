from tools.nc_reader import nc_reader
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
from utils import *
import argparse
import os
#from numpy.polynomial import Polynomial as poly

parser = argparse.ArgumentParser(description='Create power spectrum plot.')

parser.add_argument('--path',
                    type=str,
                    help='file path')

parser.add_argument('--ncfile',
                    type=str,
                    nargs='+',
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

for i, nc in enumerate(ncfile):
    ncfile[i] = os.path.join(path, nc)
    
efile = os.path.join(path, efile)

sp1 = os.path.join(path, 'spectrum_exp1_decay.asc')
sp2 = os.path.join(path, 'spectrum_exp2_decay.asc')


def find_steps():
    ts, kes, _ = np.loadtxt(efile, skiprows=1, unpack=True)
    kemax = kes.max()
    e1 = np.exp(1.0)
    e2 = np.exp(2.0)
    idx1 = find_nearest(kes, kemax / e1)
    idx2 = find_nearest(kes, kemax / e2)

    ncreader = nc_reader()
    mins = [10000, 10000]
    fn = ['', '']
    steps = [-1, -1]
    for	nc in ncfile:
        ncreader.open(nc)
        t_all = ncreader.get_all('t')

        step1 = find_nearest(t_all, ts[idx1])
        step2 = find_nearest(t_all, ts[idx2])
        ncreader.close()

        diff = abs(ts[idx1] - t_all[step1])
        if mins[0] > diff:
            mins[0] = diff
            steps[0] = step1
            fn[0] = nc
            print("Closest time to", ts[idx1], "is", t_all[steps[0]])

        diff = abs(ts[idx2] - t_all[step2])
        if mins[1] > diff:
            mins[1] = diff
            steps[1] = step2
            fn[1] = nc
            print("Closest time to", ts[idx2], "is", t_all[steps[1]])

    print("Step of exp(-1) decay:", steps[0])
    print("Step of exp(-2) decay:", steps[1])
    return steps[0], steps[1], fn[0], fn[1]

step1, step2, fn1, fn2 = find_steps()

if not os.path.exists(sp1):
    pathname = os.path.splitext(fn1)[0]
    print('genspec --filename ' + fn1 + ' --step ' + str(step1))
    os.system('genspec --filename ' + fn1 + ' --step ' + str(step1))
    os.replace(pathname + '_spectrum.asc', sp1)

if not os.path.exists(sp2):
    pathname = os.path.splitext(fn2)[0]
    print('genspec --filename ' + fn2 + ' --step ' + str(step2))
    os.system('genspec --filename ' + fn2 + ' --step ' + str(step2))
    os.replace(pathname + '_spectrum.asc', sp2)

def plot_spectrum(ax, ff, label, fit=False):
    k, p = np.loadtxt(ff, skiprows=3 ,unpack=True)

    if fit:
        hi = find_nearest(k, 300)
        plt.loglog(k[1:hi], max(p) * k[1:hi] ** (-5.0 / 3.0),
                   linestyle='dashed', color='black',
                   label=r'$P(|\bm{K}|)\propto P_{\max}|(\bm{K})|^{-5/3}$', zorder=10)

#        lo = find_nearest(k, 1)
#        hi = find_nearest(k, 250)
#        # linear fit: log10(p) = m * log10(k) + q
#        # p = 10 ** (m * log10(k) + q)
#        p_fitted = poly.fit(x=np.log10(k[lo:hi]), y=np.log10(p[lo:hi]), deg=1)
#        p_fitted = p_fitted.convert()
#        np.polynomial.set_default_printstyle('ascii')
#        print("Fitted polynomial:", p_fitted)
#        q = p_fitted.coef[0]
#        m = p_fitted.coef[1]
#        ax.loglog(k[lo:hi], 10 ** (m * np.log10(k[lo:hi]) + q), linestyle='dashed',
#                  color='black', zorder=10,
#                  label=r'$\log_{10}\mathcal{K}=' + str(round(m, 3)) + '\log_{10}|k|+'+str(round(q, 3)) + '$')

    
    ax.loglog(k, p, label=label)

mpl.rcParams['font.size'] = 10

plt.figure(figsize=(8, 2.5), dpi=400)
#fig, axs = plt.subplots(2, 1, figsize=(8, 5), dpi=400, sharex=True, sharey=False)
#grid = axs.flatten()

ax = plt.gca()

plot_spectrum(ax, sp1, label=r'$\mathcal{K}_{\max}/e$', fit=True)
plot_spectrum(ax, sp2, label=r'$\mathcal{K}_{\max}/e^2$')

ax.grid(which='both', zorder=-1)
ax.set_xlabel(r'wavenumber magnitude, $|\bm{K}| = |(\bm{k}, m)|$')
ax.set_ylabel(r'power spectrum, $P(|\bm{K}|)$')

ax.set_ylim([0.001, 10**8])
ax.set_xlim([1, 400])

ax.legend(loc='upper center', ncol=3, bbox_to_anchor=(0.5, 1.35))

plt.tight_layout()
save_figure(plt=plt, figpath=save_path, fignum=fignum, overwrite=overwrite)
plt.close()
