# Figure 6: LD decay of synonymous mutations across pops
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
            open(f"../analysis/parsed_data/{p}.unphased.outside_domains.decay.bp", "rb")
        )

ms = 3
lw = 1
markers = ["x", "+", "1"]
colors = Colorblind[8]

bins = data[p]["bins"]
bins = [(1, 10)] + bins[3:]
bin_mids = np.mean(bins, axis=1)

# combine first three bins, for distances [1, 10]

fig = plt.figure(5, figsize=(6.5, 5))

fig.clf()

i = 0
for continent, pops in populations.items():
    i += 1
    ax1 = plt.subplot(2, 3, i)
    ax2 = plt.subplot(2, 3, i + 3)
    label_syn = "Synonymous"
    label_mis = "Missense"
    for pop in pops:
        first_bins_syn = np.sum(
            list(data[pop]["ld"]["synonymous"].values())[:3], axis=0
        )
        first_bins_mis = np.sum(list(data[pop]["ld"]["missense"].values())[:3], axis=0)
        syn2 = [first_bins_syn[0] / first_bins_syn[2]] + [
            v[0] / v[2] for v in list(data[pop]["ld"]["synonymous"].values())[3:]
        ]
        mis2 = [first_bins_mis[0] / first_bins_mis[2]] + [
            v[0] / v[2] for v in list(data[pop]["ld"]["missense"].values())[3:]
        ]
        syn1 = [first_bins_syn[3] / first_bins_syn[2]] + [
            v[3] / v[2] for v in list(data[pop]["ld"]["synonymous"].values())[3:]
        ]
        mis1 = [first_bins_mis[3] / first_bins_mis[2]] + [
            v[3] / v[2] for v in list(data[pop]["ld"]["missense"].values())[3:]
        ]
        ax1.plot(
            bin_mids, syn2, "o--", color=colors[0], ms=ms, lw=lw, label=label_syn
        )
        ax1.plot(bin_mids, mis2, "s--", color=colors[1], ms=ms, lw=lw, label=label_mis)
        ax2.plot(
            bin_mids, syn1, "o--", color=colors[0], ms=ms, lw=lw, label=label_syn
        )
        ax2.plot(bin_mids, mis1, "s--", color=colors[1], ms=ms, lw=lw, label=label_mis)
        label_syn = None
        label_mis = None
    #rhos = bin_mids * 1e4 * 1e-8
    #ax1.plot(bin_mids, (5 + rhos / 2) / (11 + 13 * rhos / 2 + rhos ** 2 / 2), "k--", lw=lw)
    ax2.plot(bin_mids, 0 * bin_mids, "k--", lw=lw)
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax2.set_xscale("log")
    ax1.set_xlim(4, 5e5)
    ax2.set_xlim(4, 5e5)
    ax2.set_xlabel("Distance (bp)")
    ax1.set_ylim(1e-3, 1)
    ax2.set_ylim(-1, 1.5)
    if i == 1:
        ax1.set_ylabel("$\sigma_d^2$")
        ax2.set_ylabel("$\sigma_d^1$")
        ax1.legend()
    ax1.set_title(continent)

fig.tight_layout()

plt.savefig(f"supp_ld_decay.pdf")
