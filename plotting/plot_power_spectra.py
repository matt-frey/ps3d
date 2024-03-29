from nc_reader import nc_reader
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
from utils import *
import argparse
import os
from mpl_toolkits.axes_grid1 import ImageGrid
from numpy.polynomial import Polynomial as poly

mpl.rcParams['font.size'] = 13

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

def find_steps():
    ts, kes, _ = np.loadtxt(efile, skiprows=1, unpack=True)
    kemax = kes[0]
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

def plot_spectrum(ax, ff, label, fit=False, kbegin=None, kend=None, k53=None, surface=False):
    k, p = np.loadtxt(ff, skiprows=3 ,unpack=True)

    if fit:
        if surface:
            lab = r'$\propto|\bm{k}|^{-5/3}$'
        else:
            lab = r'$\propto|\bm{K}|^{-5/3}$'
        hi = find_nearest(k, k53)
        ax.loglog(k[1:hi], max(p) * k[1:hi] ** (-5.0 / 3.0),
                  linestyle='dashdot', color='black',
                  label=lab, zorder=10)
                  #label=r'$\propto P_{\max}|\bm{K}|^{-5/3}$', zorder=10)

        #lo = find_nearest(k, kbegin)
        #hi = find_nearest(k, kend)
        #
        #plo = p[lo]
        ## linear fit: log10(p) = log10(p[lo]) + m * log10(K)
        ## --> p = p[lo] * K^m
        #p_fitted = poly.fit(x=np.log10(k[lo:hi]), y=np.log10(p[lo:hi]), deg=1)
        #p_fitted = p_fitted.convert()
        #np.polynomial.set_default_printstyle('ascii')
        #print("Fitted polynomial:", p_fitted)
        #q = p_fitted.coef[0]
        #m = p_fitted.coef[1]
        #print("val", 10 ** q * k[lo] ** m)
        #
        #if surface:
        #    lab = r'$\propto|\bm{k}|$'
        #else:
        #    lab = r'$\propto|\bm{K}|$'
        #
        #ax.loglog(k[lo:hi], 10 ** q * k[lo:hi] ** m, linestyle='dashed',
        #          color='gray', zorder=10,
        #          label=lab + r'$^{'+str(round(m, 2)) + '}$')
#       #           label=r'$\propto' + str(round(10 ** q, 2)) + r'|\bm{K}|^{'+str(round(m, 2)) + '}$')


    ax.loglog(k, p, label=label)
    ax.legend(loc='lower left', ncol=1) #, bbox_to_anchor=(0.5, 1.35))


step1, step2, fn1, fn2 = find_steps()


lo1 = os.path.join(path, 'lower_spectrum_exp1_decay.asc')
up1 = os.path.join(path, 'upper_spectrum_exp1_decay.asc')
lo2 = os.path.join(path, 'lower_spectrum_exp2_decay.asc')
up2 = os.path.join(path, 'upper_spectrum_exp2_decay.asc')

if not os.path.exists(lo1) or not os.path.exists(up1):
    pathname = os.path.splitext(fn1)[0]
    print('genspec2d --filename ' + fn1 + ' --step ' + str(step1))
    os.system('genspec2d --filename ' + fn1 + ' --step ' + str(step1))
    os.replace(pathname + '_lower_boundary_spectrum.asc', lo1)
    os.replace(pathname + '_upper_boundary_spectrum.asc', up1)

if not os.path.exists(lo2) or not os.path.exists(up2):
    pathname = os.path.splitext(fn2)[0]
    print('genspec2d --filename ' + fn2 + ' --step ' + str(step2))
    os.system('genspec2d --filename ' + fn2 + ' --step ' + str(step2))
    os.replace(pathname + '_lower_boundary_spectrum.asc', lo2)
    os.replace(pathname + '_upper_boundary_spectrum.asc', up2)

sp1 = os.path.join(path, 'spectrum_exp1_decay.asc')
sp2 = os.path.join(path, 'spectrum_exp2_decay.asc')

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

fig = plt.figure(figsize=(13.5, 4), dpi=200)
grid = ImageGrid(fig, 111,
                nrows_ncols=(1, 3),
                aspect=False,
                axes_pad=(0.4, 0.3),
                direction='row',
                share_all=True,
                cbar_location="right",
                cbar_mode='none',
                cbar_size="4%",
                cbar_pad=0.1)

plot_spectrum(grid[0], sp1, label=r'$\mathcal{K}(t)\approx\mathcal{K}(0)/e$',
              fit=True, kbegin=2, kend=200, k53=300)
plot_spectrum(grid[0], sp2, label=r'$\mathcal{K}(t)\approx\mathcal{K}(0)/e^2$',
              fit=False, kbegin=2, kend=200, k53=300)


plot_spectrum(grid[1], lo1, label=r'$\mathcal{K}(t)\approx\mathcal{K}(0)/e$', fit=True,
              kbegin=2, kend=100, k53=100, surface=True)
plot_spectrum(grid[1], lo2, label=r'$\mathcal{K}(t)\approx\mathcal{K}(0)/e^2$', fit=False)
plot_spectrum(grid[2], up1, label=r'$\mathcal{K}(t)\approx\mathcal{K}(0)/e$', fit=True,
              kbegin=2, kend=100, k53=100, surface=True)
plot_spectrum(grid[2], up2, label=r'$\mathcal{K}(t)\approx\mathcal{K}(0)/e^2$', fit=False)

add_annotation(grid[0], r'$z = [-\pi/2, \pi/2]$', xy=(0.03, 1.06))
add_annotation(grid[1], r'$z = -\pi/2$', xy=(0.03, 1.06))
add_annotation(grid[2], r'$z =  \pi/2$', xy=(0.03, 1.06))


grid[0].set_ylabel(r'power spectrum, $P(|\bm{K}|)$ or $P(|\bm{k}|)$')

xlab = r'wavenumber magnitude, $|\bm{K}| = |(\bm{k}, m)|$'
grid[0].set_xlabel(xlab)

xlab = r'wavenumber magnitude, $|\bm{k}| = (k, l)$'
grid[1].set_xlabel(xlab)
grid[2].set_xlabel(xlab)

for i in range(3):
    grid[i].grid(which='both', zorder=-1)
    grid[i].set_ylim([1.0e-6, 0.15])
    grid[i].set_xlim([1, 300])

#plt.tight_layout()
save_figure(plt=plt, figpath=save_path, fignum=fignum, overwrite=overwrite)
plt.close()
