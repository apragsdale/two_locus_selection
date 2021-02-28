# Figure 6: LD within and between gene domains
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
            open(f"../analysis/parsed_data/{pop}.unphased.domains.bp", "rb")
        )
        data[pop] = ld_pop
        bs_pop = pickle.load(
            open(
                f"../analysis/parsed_data/{pop}.unphased.domains.sigmad1.bootstrap_std_err.bp",
                "rb",
            )
        )
        bs[pop] = bs_pop

ms = 2
lw = 0.5
markers = ["x", "+", "1"]
colors = Colorblind[8]

def plot_inside_outside(pop, ax, legend=False, ylabel=False, xlabel=True):
    # plot missense and synonymous for each, in the order:
    # gene-wide, within domains, between domains, all outside domains, matched to
    #   distances within, matched to distances between
    to_plot = []
    err = []
    # add gene-wide values
    data_all = pickle.load(
        open(f"../analysis/parsed_data/{pop}.unphased.within.between.bp", "rb")
    )
    bs_all = pickle.load(
        open(
            f"../analysis/parsed_data/{pop}.unphased.within.between.sigmad1.bootstrap_std_err.bp",
            "rb",
        )
    )
    to_plot.append(
        data_all["ld_within"]["synonymous"][3] / data_all["ld_within"]["synonymous"][2]
    )
    to_plot.append(
        data_all["ld_within"]["missense"][3] / data_all["ld_within"]["missense"][2]
    )
    err.append(bs_all["bs_within"]["synonymous"])
    err.append(bs_all["bs_within"]["missense"])
    # add within and between values
    to_plot.append(
        data[pop]["ld_within_domains"]["synonymous"][3]
        / data[pop]["ld_within_domains"]["synonymous"][2]
    )
    to_plot.append(
        data[pop]["ld_within_domains"]["missense"][3]
        / data[pop]["ld_within_domains"]["missense"][2]
    )
    to_plot.append(
        data[pop]["ld_between_domains"]["synonymous"][3]
        / data[pop]["ld_between_domains"]["synonymous"][2]
    )
    to_plot.append(
        data[pop]["ld_between_domains"]["missense"][3]
        / data[pop]["ld_between_domains"]["missense"][2]
    )
    err.append(bs[pop]["bs_within"]["synonymous"])
    err.append(bs[pop]["bs_within"]["missense"])
    err.append(bs[pop]["bs_between"]["synonymous"])
    err.append(bs[pop]["bs_between"]["missense"])
    # add all pairs outside of domains
    data_outside = pickle.load(
        open(f"../analysis/parsed_data/{pop}.unphased.outside_domains.bp", "rb")
    )
    bs_outside = pickle.load(
        open(
            f"../analysis/parsed_data/{pop}.unphased._outside_domains.sigmad1.bootstrap_std_err.bp",
            "rb",
        )
    )
    to_plot.append(
        data_outside["ld_outside_totals"]["synonymous"][3]
        / data_outside["ld_outside_totals"]["synonymous"][2]
    )
    to_plot.append(
        data_outside["ld_outside_totals"]["missense"][3]
        / data_outside["ld_outside_totals"]["missense"][2]
    )
    err.append(bs_outside["bs_outside_domains"]["synonymous"])
    err.append(bs_outside["bs_outside_domains"]["missense"])
    # add pairs outside domains matched to within domain distances
    to_plot.append(
        data_outside["ld_weighted_within"]["synonymous"][3]
        / data_outside["ld_weighted_within"]["synonymous"][2]
    )
    to_plot.append(
        data_outside["ld_weighted_within"]["missense"][3]
        / data_outside["ld_weighted_within"]["missense"][2]
    )
    err.append(bs_outside["bs_weighted_within"]["synonymous"])
    err.append(bs_outside["bs_weighted_within"]["missense"])
    # add pairs outside domains matched to between domain distances
    to_plot.append(
        data_outside["ld_weighted_between"]["synonymous"][3]
        / data_outside["ld_weighted_between"]["synonymous"][2]
    )
    to_plot.append(
        data_outside["ld_weighted_between"]["missense"][3]
        / data_outside["ld_weighted_between"]["missense"][2]
    )
    err.append(bs_outside["bs_weighted_between"]["synonymous"])
    err.append(bs_outside["bs_weighted_between"]["missense"])

    to_plot = np.array(to_plot)
    err = np.array(err)

    labels = [
        "All\nSNPs",
        "W/in\ndom.",
        "B/tw\ndom.",
        "Non-\ndom.",
        "Match\nto w/in",
        "Match\nto b/tw",
    ]
    xx_labels = np.arange(1, len(labels) + 1)
    xx = xx_labels # np.reshape(np.array([xx_labels - 0.2, xx_labels + 0.2]).T, (1, 12)).flatten()
    ax.plot((xx[0] - 1, xx[-1] + 1), (0, 0), "k--", lw=lw, label=None)

    ax.plot(
        xx - 0.2,
        to_plot.reshape((6, 2))[:, 0],
        "o",
        ms=ms,
        color=colors[0],
        label="Synonymous",
    )
    ax.plot(
        xx + 0.2,
        to_plot.reshape((6, 2))[:, 1],
        "s",
        ms=ms,
        color=colors[1],
        label="Missense",
    )

    ax.errorbar(
        xx - 0.2,
        to_plot.reshape((6, 2))[:, 0],
        yerr=1.96 * err.reshape((6, 2))[:, 0],
        linestyle="None",
        linewidth=lw,
        color="gray",
    )
    ax.errorbar(
        xx + 0.2,
        to_plot.reshape((6, 2))[:, 1],
        yerr=1.96 * err.reshape((6, 2))[:, 1],
        linestyle="None",
        linewidth=lw,
        color="gray",
    )
    ax.set_xticks(xx)
    ax.set_xticklabels(labels, fontsize=5)
    if xlabel:
        ax.set_xlabel("Mutation classes")
    if ylabel:
        ax.set_ylabel("$\sigma_d^1$")
    if legend:
        ax.legend()
    ax.set_title(pop)
    ax.set_xlim(xx[0] - 0.5, xx[-1] + 0.5)

fig = plt.figure(5, figsize=(6.5, 8))
fig.clf()
i = 0
for continent, pops in populations.items():
    for pop in pops:
        i += 1
        ax = plt.subplot(5, 3, i)
        legend = False
        if i % 3 == 1:
            ylabel = True
        else:
            ylabel = False
        if i > 12:
            xlabel = True
        else:
            xlabel = False
        plot_inside_outside(pop, ax, xlabel=xlabel, ylabel=ylabel, legend=legend)

fig.tight_layout()

plt.savefig(f"supp_domain_all_pops.pdf")
