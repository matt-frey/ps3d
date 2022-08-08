import numpy as np
import matplotlib.pyplot as plt
from utils import *
import argparse

parser = argparse.ArgumentParser(description='Create vor2vel convergence plot.')
parser.add_argument('--path',
                    type=str,
                    help='output path')
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
save_path = args.save_path
overwrite = args.overwrite
fignum = args.fignum

fig, axs = plt.subplots(1, 2, figsize=(7, 3.5), dpi=200, sharex=True)

for i in [2, 5]: #range(1, 6)
    fname = os.path.join(path, 'test_vor2vel_' + str(i) + '.asc')

    if not os.path.exists(fname):
        print("Error: You need to run 'tests/test_vor2vel' first. Exiting.")
        exit()

    nz, emax, erms = np.loadtxt(fname, skiprows=0, unpack=True)

    axs[0].plot(nz, emax, label=r'$\bm{\omega}_' + str(i) + r'$')
    axs[1].plot(nz, erms)

    if i == 2 or i == 5:
        label2 = None
        label3 = None
        if i == 2:
            label2 = r'$h^{2}$'
            label3 = r'$h^{3}$'
        h = 1.0 / np.asarray(nz)
        axs[0].plot(nz, emax[0] * (h / h[0]) ** 3, label=label3,
                    linestyle='dashed', color='black', linewidth=1)
        axs[1].plot(nz, erms[0] * (h / h[0]) ** 2, label=label2,
                    linestyle='dotted', color='black', linewidth=1)

#plt.plot(nz, err[0]* (h / h[0]) ** 2, label=r'$h^{2}$', linestyle='dashed')
#plt.plot(nz, err[0]* (h / h[0]), label=r'$h$', linestyle='dashed')


axs[0].set_xlabel(r'$n_x = n_y = n_z$')
axs[0].set_xscale('log', base=2)
axs[0].set_yscale('log', base=10)
axs[0].set_ylabel(r'$|\bm{u}_{\mathrm{ref}} - \bm{u}|_{\max}$')
axs[1].set_xlabel(r'$n_x = n_y = n_z$')
axs[1].set_xscale('log', base=2)
axs[1].set_yscale('log', base=10)
axs[1].set_ylabel(r'$|\bm{u}_{\mathrm{ref}} - \bm{u}|_{\mathrm{rms}}$')
plt.figlegend(loc='upper center', ncol=5, bbox_to_anchor=(0.5, 1.0))
plt.tight_layout()
fig.subplots_adjust(top=0.88)

save_figure(plt=plt, figpath=save_path, fignum=fignum, overwrite=overwrite)
plt.close()
