import os
import numpy as np
import matplotlib.pyplot as plt
from collections import OrderedDict
from plot_helpers import ls, moving_average

labels = {
        'base_elasticity_cr_False': r'$\mathring{H}(\mathrm{sym})$',
        'base_elasticity_cr_True': r'$CR(10^{-2})+\mathring{H}(\mathrm{sym})$',
        'base_laplace_cr_False': r'$\mathring{H}^1$',
        'base_laplace_cr_True': r'$CR(10^{-2})+\mathring{H}^1$',
        }
m = [(0.05, 0.1), (0.05, 0.1), (0.00, 0.1), (0.00, 0.1)]
lts = [ls["densely dotted"], ls["densely dotted"], ls["densely dashed"], ls["densely dashed"]]
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']
colors = [colors[0], colors[1], colors[2], colors[3]]

plt.rcParams['text.latex.preamble'] = [r"\usepackage{lmodern}"]
params = {'text.usetex': True,
          'font.size': 7,
          'font.family': 'lmodern',
          'text.latex.unicode': True,
          'axes.titlesize': 11,
          }
plt.rcParams.update(params)
# plt.rc('text', usetex=True)
# plt.rc('font', family='serif')
f, axarr = plt.subplots(1, 2, figsize=(15.24/2.54, 5.6/2.54))
for counter, name in enumerate([ 'base_elasticity_cr_False','base_elasticity_cr_True',
              'base_laplace_cr_False', 'base_laplace_cr_True']):
    d = "./output/stokes/" + name + "/"
    Jvals = np.load(d + "Jvals.npy")
    cnorms = np.load(d + "cnorms.npy")
    gnorms = np.load(d + "gnorms.npy")
    pde_solves = np.load(d + "pde_solves.npy")
    axarr[0].plot(pde_solves, Jvals, label=labels[name], linewidth=1.0, color=colors[counter], linestyle=lts[counter])
    # axarr[1].semilogy(pde_solves, gnorms/gnorms[0], label=labels[name], marker=ms[counter],
    #                   markevery=m[counter], markersize=2, linewidth=0.5)
    axarr[1].semilogy(pde_solves, gnorms, label=labels[name], linewidth=1.0, color=colors[counter], linestyle=lts[counter])
    # axarr[1].semilogy(moving_average(pde_solves, n=10), moving_average(cnorms, n=10), label=labels[name], marker=ms[counter],
                      # markevery=m[counter], markersize=2, linewidth=0.5)
    # axarr[0].semilogy(pde_solves, cnorms, label=labels[name])

for ax in axarr:
    ax.set_xlabel("PDE Solves")

axarr[0].legend()
axarr[0].set_title('Functional value')
axarr[1].legend()
axarr[1].set_title('Norm of gradient of Lagrangian')
# axarr[2].legend()
# axarr[2].set_title('Constraint violation')
mydir = r'./output/stokes/img/'
plt.tight_layout(pad=0.)
f.savefig(mydir + "stokes_convergence.pdf", dpi=1000)
