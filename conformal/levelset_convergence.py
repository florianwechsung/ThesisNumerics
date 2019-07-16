import numpy as np
import matplotlib
matplotlib.use('Agg')  # noqa
import matplotlib.pyplot as plt
import os
from plot_helpers import ls, moving_average


plt.figure()
outdirs = [
    "levelset-base-elasticity-cr-None",
    "levelset-base-elasticity-cr-0.3",
    "levelset-base-elasticity-cr-0.01",
    "levelset-base-elasticity-cr-0.001",
    "levelset-base-laplace-cr-0.3",
    "levelset-base-laplace-cr-0.01",
    "levelset-base-laplace-cr-0.001",
]
labels = [
    r"$\mathring{H}^1(\mathrm{sym})$",
    r"$CR(0.3)+\mathring{H}^1(\mathrm{sym})$",
    r"$CR(10^{-2})+\mathring{H}^1(\mathrm{sym})$",
    r"$CR(10^{-3})+\mathring{H}^1(\mathrm{sym})$",
    r"$CR(0.3)+{H}^1$",
    r"$CR(10^{-2})+\mathring{H}^1$",
    r"$CR(10^{-3})+\mathring{H}^1$"
]

plt.figure(figsize=(15.24/2.54, 7.6/2.54))
plt.rcParams['text.latex.preamble'] = [r"\usepackage{lmodern}"]
params = {'text.usetex': True,
          'font.size': 7,
          'font.family': 'lmodern',
          'text.latex.unicode': True,
          'axes.titlesize': 11,
          }
plt.rcParams.update(params)
ax = plt.subplot(111)
ax.set_title("Relative gradient norm")
ax.set_xlabel("Iteration")
lts = [ls["densely dotted"], ls["densely dotted"], ls["densely dotted"], ls["densely dotted"], ls["densely dashed"], ls["densely dashed"], ls["densely dashed"]]
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']
colors = [colors[0], colors[1], colors[2], colors[3], colors[1], colors[2], colors[3]]

tablecolumns = []
header = []
for i, outdir in enumerate(outdirs):
    if os.path.exists(os.path.join("./output/levelset", outdir)):
        fi = "gradient_norms"
        # fi = "boundary_derivatives"
        # fi = "objective_values"
        y = np.loadtxt(os.path.join("output", "levelset", outdir, fi + ".txt"))
        y = moving_average(y, n=5)
        if fi == "objective_values":
            shift = 20
            y = y[shift:]
            x = list(range(20, len(y)+20))
            plt.plot(x, y, label=labels[i], linestyle=lts[i])
        else:
            y = y[1:]/y[1]
            x = list(range(0, len(y)))
            plt.semilogy(x, y, label=labels[i], linestyle=lts[i], color=colors[i])
        tablecolumns.append(x)
        tablecolumns.append(y)
        header.append(outdir + "-x")
        header.append(outdir + "-y")

plotdir = "output/levelset/img/"
os.makedirs(plotdir, exist_ok=True)
plt.legend(loc='lower left')
plt.tight_layout(pad=0.)
plt.savefig(os.path.join(plotdir, "levelset_convergence.pdf"))

np.savetxt(plotdir + "levelsetconvergencedata.txt", np.vstack(tablecolumns).T, delimiter=",", newline="\n", header=",".join(header), comments='')
