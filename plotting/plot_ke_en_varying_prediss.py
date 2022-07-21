import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams.update({
    "figure.figsize": (9, 6),
    "figure.dpi": 200,
    "font.family": "serif",
    "font.size": 20,
    "text.usetex": True,
    'text.latex.preamble': "\n".join([
        r"\usepackage{amsmath}",
        r"\usepackage[utf8]{inputenc}",
        r"\usepackage[T1]{fontenc}",
#        r"\usepackage{siunitx}",
        ])
})


labels = [r'$C = 10$',
          r'$C = 20$',
          r'$C = 30$',
          r'$C = 60$',
          r'$C = 100$'
]

branches = ['pred10',
            'pred20',
            'pred30',
            'pred60',
            'pred100'
]

grids = [32]

mpl.rcParams['font.size'] = 11

fig, axs = plt.subplots(2, 1, figsize=(7, 2*2.5), dpi=200, sharex=True)

i = 0

for branch in branches:
        for grid in grids:
                t, ke, en = np.loadtxt('paper_runs/beltrami_' + str(grid) + '_' + branch + '_ecomp.asc',
                                       skiprows=1, unpack=True)

                ncelli = 1.0 / grid ** 3

                # calculate mean KE and mean EN
                ke *= ncelli
                en *= ncelli
                
                print("initial <KE>", ke[0])
                print("initial <EN>", en[0])

                if labels[i]:
                        label = labels[i]
                        i = i + 1
                else:
                        label = branch
                
                axs[0].plot(t, ke, label=label)
                
                axs[1].plot(t, en, label=label)

axs[1].set_xlabel(r'time, $t$')
axs[1].set_ylabel(r'average enstrophy, $\langle\Upsilon\rangle$')

axs[0].set_ylabel(r'average kinetic energy, $\langle\mathcal{K}\rangle$')

axs[0].legend(loc='upper center', ncol=5, bbox_to_anchor=(0.5, 1.2))

plt.tight_layout()

plt.savefig('fig1.eps', format='eps') #, bbox_inches='tight', dpi=200)

plt.close()
