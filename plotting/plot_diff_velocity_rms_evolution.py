from nc_reader import nc_reader
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
from utils import *
import argparse
import os
from scipy.optimize import curve_fit

parser = argparse.ArgumentParser(description='Create cross section figures.')
parser.add_argument('--filename',
                    type=str,
                    help='output file')

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
fname = args.filename
save_path = args.save_path
overwrite = args.overwrite
fignum = args.fignum

print()
print("\tFilename:    ", fname)
print("\tSave path:   ", save_path)
print("\tOverwrite:   ", overwrite)
print("\tFignum:      ", fignum)
print()

ncr = nc_reader()

ncr.open(fname)

t = ncr.get_all('t')
n = len(t)

u0 = ncr.get_dataset(0, 'x_velocity')
v0 = ncr.get_dataset(0, 'y_velocity')
w0 = ncr.get_dataset(0, 'z_velocity')

# exclude first entry since it is obviously zero --> vector length n-1
diff_mag_rms = np.zeros(n-1)
t = t[1:]
for step in range(1, n):
    u = ncr.get_dataset(step, 'x_velocity')
    v = ncr.get_dataset(step, 'y_velocity')
    w = ncr.get_dataset(step, 'z_velocity')

    diff_mag = np.sqrt((u - u0) ** 2 + (v - v0) ** 2 + (w - w0) ** 2)
    diff_mag_rms[step-1] = np.sqrt((diff_mag ** 2).mean(axis=(0, 1, 2)))

ncr.close()

# take natural logarithm
log_diff_mag_rms = np.log(diff_mag_rms)


idx1 = find_nearest(t, 53)
idx2 = find_nearest(t, 58)

print("Fit start: ", t[idx1])
print("Fit end:   ", t[idx2])
print("Num points:", idx2 - idx1 + 1)

def linear_fit_func(x, a, b):
    return a * x + b

xdata = t[idx1:idx2+1]
ydata = log_diff_mag_rms[idx1:idx2+1]


popt, pcov = curve_fit(linear_fit_func, xdata, ydata, method='lm', ftol=1e-8,
                       absolute_sigma=True)


ypred = linear_fit_func(xdata, *popt)

#
# Bootstrap growth rate and calculate standard deviation:
#
n_sample = 10000

np.random.seed(seed=42)

def get_sample(xd, yd):
    indices = np.random.randint(low=0, high=len(xd), size=len(xd))
    return np.take(xd, indices), np.take(yd, indices)


aas = np.zeros(n_sample)
bbs = np.zeros(n_sample)
for i in range(n_sample):
    x, y = get_sample(xdata, ydata)
    popt, pcov = curve_fit(linear_fit_func, x, y, method='lm', ftol=1e-8, absolute_sigma=True)
    aas[i], bbs[i] = popt

a_std = aas.std()
b_std = bbs.std()

a, b = popt
print("Fit: y = a x + b")
print("Coefficients a =", a, "b =", b)
print("Std. dev. of a = ", a_std, "and of b = ", b_std)

plt.figure(figsize=(7, 3), dpi=200)
plt.semilogy(xdata, np.exp(ypred), 'k--', base=np.e,
             label=r'$|\Delta\bm{u}|_{\mathrm{rms}}\propto(' + str(round(a, 4)) + \
             '\pm' + str(round(a_std, 4)) + ')t$')

plt.semilogy(t, diff_mag_rms, base=np.e, color=colors[0])
plt.grid(zorder=-2)
plt.xlabel(r'time, $t$')
plt.ylabel(r'$|\Delta\bm{u}|_{\mathrm{rms}}$')

plt.yticks([np.exp(-6), np.exp(-5), np.exp(-4), np.exp(-3), np.exp(-2), np.exp(-1)],
           [r'$e^{-6}$', r'$e^{-5}$', r'$e^{-4}$', r'$e^{-3}$', r'$e^{-2}$', r'$e^{-1}$'])

plt.legend(loc='upper left')
plt.tight_layout()
save_figure(plt=plt, figpath=save_path, fignum=fignum, overwrite=overwrite)
plt.close()
