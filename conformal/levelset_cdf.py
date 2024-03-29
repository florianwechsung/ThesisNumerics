import os
import numpy as np
import matplotlib.pyplot as plt
from collections import OrderedDict
from plot_helpers import ls
labels = {
        'base-elasticity-cr-None': r'$\mathring{H}(\mathrm{sym})$',
        'base-elasticity-cr-0p3': r'$CR(0.3)+\mathring{H}(\mathrm{sym})$',
        'base-elasticity-cr-0p01': r'$CR(10^{-2})+\mathring{H}(\mathrm{sym})$',
        'base-laplace-cr-0p3': r'$CR(0.3)+\mathring{H}^1$',
        'base-laplace-cr-0p01': r'$CR(10^{-2})+\mathring{H}^1$',
        }
ms = ['o', 'o', 'o', 's', 's']
m = [(0.00, 0.1), (0.00, 0.1), (0.00, 0.1), (0.05, 0.1), (0.05, 0.1), (0.00, 0.1)]
lts = [ls["densely dotted"], ls["densely dotted"], ls["densely dotted"], ls["densely dashed"], ls["densely dashed"]]
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']
colors = [colors[0], colors[1], colors[2], colors[1], colors[2], colors[4]]
plt.figure(figsize=(15.24/2.54, 7.6/2.54))
plt.rcParams['text.latex.preamble'] = [r"\usepackage{lmodern}"]
params = {'text.usetex': True,
          'font.size': 7,
          'font.family': 'lmodern',
          'text.latex.unicode': True,
          'axes.titlesize': 11,
          }
plt.rcParams.update(params)
mydir = r'./output/levelset/img/'
counter = 0
qs = [1.5, 1.75, 2.0, 2.5, 3.0]
print("&" + '&'.join(str(q) for q in qs) + r"&\text{max}\\\hline")

tablecolumns = []
header = []
for key in labels.keys():
    f = "inv_cdf_levelset-" + key + ".npy"
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
    if counter == 5:
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
plt.savefig(mydir + "levelset_cdf.pdf", dpi=300)
qualis = np.linspace(1.0, 3.0, num=100)
for i in range(len(tablecolumns)//2):
    tablecolumns[2*i+1] = np.interp(qualis, tablecolumns[2*i], tablecolumns[2*i+1])
    tablecolumns[2*i] = qualis
np.savetxt(mydir + "levelsetcdfdata.txt", np.vstack(tablecolumns).T, delimiter=",", newline="\n", header=",".join(header), comments='')
