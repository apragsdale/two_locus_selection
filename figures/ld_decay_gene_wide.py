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

populations = {
    "Africa": ["ESN", "GWD", "LWK", "MSL", "YRI"],
    "Europe": ["CEU", "FIN", "GBR", "IBS", "TSI"],
    "East Asia": ["CDX", "CHB", "CHS", "JPT", "KHV"],
}

for c, ps in populations.items():
    for p in ps:
        data[p] = pickle.load(
            open(f"../analysis/parsed_data_v2/{p}.unphased.all.bp", "rb")
        )

ms = 2
lw = 0.5
markers = ["x", "+", "1"]
colors = Colorblind[8]

bins = data[p]["bin_edges"]
bin_mids = np.mean(bins, axis=1)

fig = plt.figure(5, figsize=(6.5, 3.5))

fig.clf()

i = 0
for continent, pops in populations.items():
    i += 1
    ax1 = plt.subplot(2, 3, i)
    ax2 = plt.subplot(2, 3, i + 3)
    label_syn = "Synonymous"
    label_mis = "Missense"
    ax2.plot(bin_mids[:-1], 0 * bin_mids[:-1], "k--", lw=lw)
    for pop in pops:
        sd2_syn = data[pop]["bins"]["all_sites"]["synonymous"]["sd2"]
        sd2_mis = data[pop]["bins"]["all_sites"]["missense"]["sd2"]
        sd1_syn = data[pop]["bins"]["all_sites"]["synonymous"]["sd1"]
        sd1_mis = data[pop]["bins"]["all_sites"]["missense"]["sd1"]
        ax1.plot(
                bin_mids[:-1], sd2_syn[:-1], "o--", color=colors[0], ms=ms, lw=lw, label=label_syn
        )
        ax1.plot(bin_mids[:-1], sd2_mis[:-1], "s--", color=colors[1], ms=ms, lw=lw, label=label_mis)
        ax2.plot(
                bin_mids[:-1], sd1_syn[:-1], "o--", color=colors[0], ms=ms, lw=lw, label=label_syn
        )
        ax2.plot(bin_mids[:-1], sd1_mis[:-1], "s--", color=colors[1], ms=ms, lw=lw, label=label_mis)
        label_syn = None
        label_mis = None
    #rhos = bin_mids * 1e4 * 1e-8
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax2.set_xscale("log")
    #ax1.set_xlim(4, 5e5)
    #ax2.set_xlim(4, 5e5)
    ax2.set_xlabel("Distance (bp)")
    #ax1.set_ylim(1e-3, 1)
    #ax2.set_ylim(-1, 1.5)
    if i == 1:
        ax1.set_ylabel("$\sigma_d^2$")
        ax2.set_ylabel("$\sigma_d^1$")
        ax1.legend()
    ax1.set_title(continent)

fig.tight_layout()

plt.savefig(f"ld_decay_gene_wide.pdf")
