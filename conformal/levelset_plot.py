import numpy as np
import matplotlib
matplotlib.use('Agg')  # noqa
import matplotlib.pyplot as plt
import os
from plot_helpers import ls


plt.figure()
outdirs = [
    "levelset-base-elasticity-cr-None",
    "levelset-base-elasticity-cr-0.3",
    "levelset-base-elasticity-cr-0.001",
    "levelset-base-laplace-cr-0.3",
    "levelset-base-laplace-cr-0.001",
]
labels = [
    r"$\mathring{H}^1(\mathrm{sym})$",
    r"$CR(0.3)+\mathring{H}^1(\mathrm{sym})$",
    r"$CR(10^{-3})+\mathring{H}^1(\mathrm{sym})$",
    r"$CR(0.3)+{H}^1$",
    r"$CR(10^{-3})+\mathring{H}^1$"
]

fig = plt.figure(figsize=(10, 4), dpi=300)
ax = plt.subplot(111)
ax.set_title("Gradient norm")
lts = [ls["densely dashed"], ls["densely dotted"], ls["dashed"], ls["densely dashdotted"], ls["solid"]]

for i, outdir in enumerate(outdirs):
    if os.path.exists(outdir):
        fi = "gradient_norms"
        # fi = "boundary_derivatives"
        # fi = "objective_values"
        y = np.loadtxt(os.path.join(outdir, fi + ".txt"))
        if fi == "objective_values":
            shift = 20
            y = y[shift:]
            x = list(range(20, len(y)+20))
            plt.plot(x, y, label=labels[i], linestyle=lts[i])
        else:
            y = y[1:]/y[1]
            x = list(range(0, len(y)))
            plt.semilogy(x, y, label=labels[i], linestyle=lts[i])

plotdir = "levelset"
os.makedirs(plotdir, exist_ok=True)
plt.legend()
plt.savefig(os.path.join(plotdir, fi + ".pdf"), bbox_inches='tight')
