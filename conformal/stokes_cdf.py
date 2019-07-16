import os
import numpy as np
import matplotlib.pyplot as plt
from plot_helpers import ls
labels = {
        'base_elasticity_cr_False': r'$\mathring{H}(\mathrm{sym})$',
        'base_elasticity_cr_True': r'$CR(10^{-2})+\mathring{H}(\mathrm{sym})$',
        'base_laplace_cr_False': r'$\mathring{H}^1$',
        'base_laplace_cr_True': r'$CR(10^{-2})+\mathring{H}^1$',
        }

m = [(0.05, 0.1), (0.05, 0.1), (0.00, 0.1), (0.00, 0.1)]
ms = ["s", "s", "o", "o"]
lts = [ls["densely dotted"], ls["densely dotted"], ls["densely dashed"], ls["densely dashed"]]
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']
colors = [colors[0], colors[1], colors[2], colors[3], colors[4]]

plt.figure(figsize=(15.24/2.54, 7.6/2.54))
plt.rcParams['text.latex.preamble'] = [r"\usepackage{lmodern}"]
params = {'text.usetex': True,
          'font.size': 7,
          'font.family': 'lmodern',
          'text.latex.unicode': True,
          'axes.titlesize': 11,
          }
plt.rcParams.update(params)
mydir = r'./output/stokes/img/'
counter = 0
qs = [1.5, 1.75, 2.0, 2.5, 3.0]
print("&" + '&'.join(str(q) for q in qs) + r"&\text{max}\\\hline")

tablecolumns = []
header = []
for key in labels.keys():
    f = "inv_cdf_" + key + ".npy"
    label = labels[key]
    xy = np.load(mydir + f)
    qual = xy[:, 1]
    x = qual
    y = xy[:, 2]
    plt.plot(x, y, label=label, marker=ms[counter],
             markevery=m[counter], markersize=4, color=colors[counter], linewidth=1, linestyle=lts[counter])
    print(label + "&"
          + "&".join(f"{100.*np.sum(qual>q)/len(qual):.1f}\%" for q in qs)
          + "&" + f"{np.max(qual):.1f}" + r"\\")
    counter += 1
    tablecolumns.append(x)
    tablecolumns.append(y)
    header.append(key + "-x")
    header.append(key + "-y")
    if counter == 4:
        qual = xy[:, 0]
        plt.plot(qual, xy[:, 2], label="Initial mesh", color=colors[counter])
        print("Initial mesh" + "&" +
              "&".join(f"{100.*np.sum(qual>q)/len(qual):.1f}\%" for q in qs)
              + "&" + f"{np.max(qual):.1f}" + r"\\")
        counter += 1
        tablecolumns.append(qual)
        tablecolumns.append(xy[:, 2])
        header.append("initial" + "-x")
        header.append("initial" + "-y")
plt.xlim((1, 1.5))
plt.xlabel(r"$\eta$")
plt.ylabel(r"Fraction of cells with $\eta(K)\le \eta$")
plt.legend()
plt.title("Mesh quality CDF")
plt.tight_layout(pad=0.)
plt.savefig(mydir + "stokes_cdf.pdf", dpi=300)
qualis = np.linspace(1.0, 3.0, num=100)
for i in range(len(tablecolumns)//2):
    tablecolumns[2*i+1] = np.interp(qualis, tablecolumns[2*i], tablecolumns[2*i+1])
    tablecolumns[2*i] = qualis
np.savetxt(mydir + "stokescdfdata.txt", np.vstack(tablecolumns).T, delimiter=",", newline="\n", header=",".join(header), comments='')
