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

parser.add_argument('--labels',
                    type=str,
                    nargs='+',
                    help='legend labels')
args = parser.parse_args()
path = args.path
labels = args.labels
save_path = args.save_path
overwrite = args.overwrite
fignum = args.fignum

fig, axs = plt.subplots(1, 2, figsize=(7, 3.5), dpi=200, sharey=True)

for i in range(len(labels)):
    fname = os.path.join(path, 'q_profile_filter_' + str(i+1) + '.asc')

    if not os.path.exists(fname):
        print("Error: You need to run 'tests/test_filter' first. Exiting.")
        exit()

    z, qs, q, diff = np.loadtxt(fname, skiprows=1, unpack=True)

    axs[0].plot(diff, z, label=labels[i])

    fname = os.path.join(path, 'q_profile_hyper_' + str(i+1) + '.asc')
    if not os.path.exists(fname):
        print("Error: You need to run 'tests/test_hyper' first. Exiting.")
        exit()
        
    z, qs, q, diff = np.loadtxt(fname, skiprows=1, unpack=True)
    axs[1].plot(diff, z)


remove_yticks(axs[1])

axs[0].ticklabel_format(style='sci', axis='x', scilimits=(0,0))
axs[1].ticklabel_format(style='sci', axis='x', scilimits=(0,0))

add_annotation(axs[0], r'filtering', xy=(0.03, 1.06), fontsize=12)
add_annotation(axs[1], r'hyperdiffusion', xy=(0.62, 1.06), fontsize=12)

axs[0].set_yticks(np.pi*np.array([-1/2, -1/4, 0, 1/4, 1/2]),
                  [r'$-\pi/2$',r'$-\pi/4$',r'$0$',r'$\pi/4$',r'$\pi/2$'])

axs[0].grid(which='both')
axs[1].grid(which='both')
axs[0].set_xlabel(r'$q_f - q_0$')
axs[0].set_ylabel(r'$z$')
axs[1].set_xlabel(r'$q_h - q_0$')
plt.figlegend(loc='upper center', ncol=5, bbox_to_anchor=(0.55, 1.0))
plt.tight_layout()
fig.subplots_adjust(top=0.88)

save_figure(plt=plt, figpath=save_path, fignum=fignum, overwrite=overwrite)
plt.close()
