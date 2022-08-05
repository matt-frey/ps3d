import numpy as np
import matplotlib.pyplot as plt
from utils import *
import argparse
import os

parser = argparse.ArgumentParser(
    description='Plot kinetic energy and enstrophy for different hyperdiffusion prefactors.')

parser.add_argument('--filepath',
                    type=str,
                    help='path to files')

parser.add_argument('--save_path',
                    type=str,
                    help='where to save the figures',
                    default=os.getcwd())

parser.add_argument('--overwrite',
                    help='overwrite figures',
                    action='store_true')

parser.add_argument('--fignum',
                    type=int,
                    default=2,
                    help='figure number')

args = parser.parse_args()
fpath = args.filepath
save_path = args.save_path
overwrite = args.overwrite
fignum = args.fignum

print()
print("\tFilepath:  ", fpath)
print("\tSave path: ", save_path)
print("\tOverwrite: ", overwrite)
print("\tFignum:    ", fignum)
print()


labels = [r'$C = 10$',
          r'$C = 20$',
          r'$C = 30$',
          r'$C = 60$',
          r'$C = 100$'
]

prediss = ['pred10',
           'pred20',
           'pred30',
           'pred60',
           'pred100'
]

grid = 32

fig, axs = plt.subplots(2, 1, figsize=(7, 4), dpi=200, sharex=True)

i = 0

for pred in prediss:
    t, ke, en = np.loadtxt(os.path.join(fpath, 'beltrami_' + str(grid) + '_' + pred + '_ecomp.asc'),
                           skiprows=1, unpack=True)

    voli = 1.0 / np.pi ** 3
    ncelli = 1.0 / grid ** 3

    # calculate mean KE and mean EN
    ke *= voli * ncelli
    en *= voli * ncelli

    #print("initial <KE>", ke[0])
    #print("initial <EN>", en[0])

    if labels[i]:
        label = labels[i]
        i = i + 1
    else:
        label = pred

    axs[0].plot(t, ke, label=label)

    axs[1].plot(t, en, label=label)

axs[1].set_xlabel(r'time, $t$')
axs[1].set_ylabel(r'enstrophy, $\Upsilon$')

axs[0].tick_params(axis='x', which='both', length=0)

axs[0].grid(zorder=-1)
axs[1].grid(zorder=-1)

axs[0].set_xlim([-1, 101])
axs[1].set_xlim([-1, 101])

axs[1].set_yticks([0, 2.5, 5, 7.5])

axs[0].set_ylabel(r'kinetic energy, $\mathcal{K}$')

axs[0].legend(loc='upper center', ncol=5, bbox_to_anchor=(0.5, 1.3))

plt.tight_layout()

save_figure(plt=plt, figpath=save_path, fignum=fignum, overwrite=overwrite)
plt.close()
