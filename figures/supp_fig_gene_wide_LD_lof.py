# Figure 5: LD within and between genes
import sys
import moments.TwoLocus
import matplotlib.pylab as plt
import numpy as np
import pickle

# set font sizes
import matplotlib

plt.rcParams["legend.title_fontsize"] = "xx-small"
matplotlib.rc("xtick", labelsize=6)
matplotlib.rc("ytick", labelsize=6)
matplotlib.rc("axes", labelsize=7)
matplotlib.rc("axes", titlesize=7)
matplotlib.rc("legend", fontsize=6)

from bokeh.palettes import Colorblind

data = {}

bs = {}

populations = {
    "Africa": ["ESN", "GWD", "LWK", "MSL", "YRI"],
    "Europe": ["CEU", "FIN", "GBR", "IBS", "TSI"],
    "East Asia": ["CDX", "CHB", "CHS", "JPT", "KHV"],
}

for continent, pops in populations.items():
    for pop in pops:
        ld_pop = pickle.load(
            open(f"../analysis/parsed_data/{pop}.unphased.within.between.bp", "rb")
        )
        data[pop] = ld_pop
        bs_pop = pickle.load(
            open(
                f"../analysis/parsed_data/{pop}.unphased.within.between.sigmad1.bootstrap_std_err.bp",
                "rb",
            )
        )
        bs[pop] = bs_pop

ms = 2
lw = 0.5

markers = ["x", "+", "1"]
colors = Colorblind[8]

## plots between genes

def plot_between(pop, data, bs, ax, legend=False, ylabel=False, xlabel=False):
    annot = "loss_of_function"
    marker = "X"
    color = colors[3]
    to_plot = [data[pop]["ld_within"][annot][3] / data[pop]["ld_within"][annot][2]]
    bins = sorted(data[pop]["ld_between"][annot].keys())
    to_plot = to_plot + [
        data[pop]["ld_between"][annot][b][3] / data[pop]["ld_between"][annot][b][2]
        for b in bins
    ]
    err = [1.96 * bs[pop]["bs_within"][annot]] + [
        1.96 * bs[pop]["bs_between"][annot][b] for b in bins
    ]
    xx = np.arange(len(to_plot))
    labels = ["w/in gene"] + [
        f"{int(b[0] / 1000)}-{int(b[1] / 1000)}kb" for b in bins
    ]
    ax.plot(
        xx, to_plot, marker + "-", ms=ms, color=color, lw=lw, label="Loss of function"
    )
    ax.errorbar(
        xx, to_plot, yerr=err, linestyle="None", linewidth=lw, color="gray"
    )
    ax.set_xticks(xx)
    ax.set_xticklabels(labels, rotation=45, ha="right")
    if xlabel:
        ax.set_xlabel("Distance between bp")
    if legend:
        ax.legend(loc="upper right")
    ax.plot(xx, np.zeros(len(xx)), "k--", lw=lw, label=None)
    if ylabel:
        ax.set_ylabel("$\sigma_d^1$")
    ax.set_title(pop)


fig = plt.figure(51, figsize=(6.5, 8))
fig.clf()

ii = 0
for continent in populations.keys():
    for pop in populations[continent]:
        ii += 1
        ax = plt.subplot(5, 3, ii)
        plot_between(pop, data, bs, ax, legend=False, ylabel=False, xlabel=False)
        if ii % 3 == 1:
            ax.set_ylabel("$\sigma_d^1$")
        if ii > 12:
            ax.set_xlabel("Distance between bp")


fig.tight_layout()
plt.savefig(f"supp_gene_wide_all_pops_lof.pdf")
