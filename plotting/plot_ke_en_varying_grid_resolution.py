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


labels = [r'$32^3$',
          r'$64^3$',
          r'$128^3$',
          r'$256^3$'
]

grids = [32, 64, 128, 256]

mpl.rcParams['font.size'] = 11

fig, axs = plt.subplots(2, 1, figsize=(7, 2*2.5), dpi=200, sharex=True)

i = 0


for grid in grids:
        t, ke, en = np.loadtxt('paper_runs/beltrami_' + str(grid) + '_ecomp.asc',
                               skiprows=1, unpack=True)
        
        ncelli = 1.0 / grid ** 3
        
        # calculate mean KE and mean EN
        ke *= ncelli
        en *= ncelli
        
        print("initial <KE>", ke[0])
        print("initial <EN>", en[0])
        
        label = labels[i]
        i = i + 1
                
        axs[0].plot(t, ke, label=label)        
        axs[1].plot(t, en, label=label)

axs[1].set_xlabel(r'time, $t$')
axs[1].set_ylabel(r'average enstrophy, $\langle\Upsilon\rangle$')

axs[0].set_ylabel(r'average kinetic energy, $\langle\mathcal{K}\rangle$')

axs[0].legend(loc='upper center', ncol=5, bbox_to_anchor=(0.5, 1.2))

plt.tight_layout()

plt.savefig('fig2.eps', format='eps') #, bbox_inches='tight', dpi=200)

plt.close()
