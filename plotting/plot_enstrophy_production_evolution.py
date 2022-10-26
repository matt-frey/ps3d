from nc_reader import nc_reader
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
from utils import *
import argparse
import os

parser = argparse.ArgumentParser(description='Plot enstrophy production evolution.')
parser.add_argument('--filename',
                    type=str,
                    help='output file')

parser.add_argument('--restartfile',
                    type=str,
                    help="another file restarted from the first one")

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
restart_file = args.restartfile
save_path = args.save_path
overwrite = args.overwrite
fignum = args.fignum

print()
print("\tFilename:    ", fname)
print("\tRestart file:", restart_file)
print("\tSave path:   ", save_path)
print("\tOverwrite:   ", overwrite)
print("\tFignum:      ", fignum)
print()

def fill_steps(data, j, lo, hi):

    t_all = data[:, 0]

    for step in range(lo, hi):
        up1[j] = data[step, 1]
        up2[j] = data[step, 2]
        up3[j] = data[step, 3]

        t[j] = t_all[step]

        j = j + 1


data_1 = np.loadtxt(fname, skiprows=1)

t1 = data_1[:, 0]
n = len(t1)

if not restart_file is None:
    data_2 = np.loadtxt(restart_file, skiprows=1)
    t2 = data_2[:, 0]
    lo1 = find_nearest(t1, t2[0])
    hi1 = find_nearest(t1, t2[-1])
    n = len(t1[0:lo1]) + len(t2) + len(t1[hi1:])

t = np.zeros(n)
up1 = np.zeros(n)
up2 = np.zeros(n)
up3 = np.zeros(n)

if restart_file is None:
    fill_steps(data_1, 0, 0, n)
else:
    fill_steps(data_1, 0, 0, lo1)
    fill_steps(data_2, lo1, 0, len(t2))
    fill_steps(data_1, lo1+len(t2), hi1, len(t1))

fig, ax = plt.subplots(1, 1, figsize=(7, 2.5), dpi=200)

#ax.axvspan(xmin=t2[0], xmax=t2[-1], color='lightgrey', zorder=-1)
ax.plot(t, up1, color=colors[0], label=r'$\upsilon_1$')
ax.plot(t, up2, color=colors[1], label=r'$\upsilon_2$')
ax.plot(t, up3, color=colors[2], label=r'$\upsilon_3$')
ax.plot(t, up1 + up2 + up3, color='black', linestyle='dashed',
        label=r'$\upsilon = \upsilon_1 + \upsilon_2 + \upsilon_3$')
ax.set_xlim([50, 70])
ax.grid(zorder=-2)

#umin = min(up1.min(), up2.min(), up3.min())
#umax = max(up1.max(), up2.max(), up3.max())
ax.set_yticks([-270, -135, 0, 135, 270, 405, 500])

ax.legend(loc='upper center', ncol=4, bbox_to_anchor=(0.5, 1.35))
ax.set_xlabel(r'time, $t$')

plt.tight_layout()
save_figure(plt=plt, figpath=save_path, fignum=fignum, overwrite=overwrite)
plt.close()
