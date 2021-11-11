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

fig = plt.figure(5, figsize=(6.5, 4.5))
fig.clf()


# plot all stats within genes
ax1 = plt.subplot2grid((5, 3), (0, 0), rowspan=3, colspan=3)

xx = np.concatenate(
    (np.linspace(1, 20, 20), np.linspace(23, 42, 20), np.linspace(45, 64, 20))
)
yy = []
err = []

for continent in ["Africa", "Europe", "East Asia"]:
    for pop in populations[continent]:
        yy.append(
            data[pop]["ld_within_domains"]["synonymous"][3]
            / data[pop]["ld_within_domains"]["synonymous"][2]
        )
        yy.append(
            data[pop]["ld_within_domains"]["missense"][3]
            / data[pop]["ld_within_domains"]["missense"][2]
        )
        yy.append(
            data[pop]["ld_between_domains"]["synonymous"][3]
            / data[pop]["ld_between_domains"]["synonymous"][2]
        )
        yy.append(
            data[pop]["ld_between_domains"]["missense"][3]
            / data[pop]["ld_between_domains"]["missense"][2]
        )
        err.append(1.96 * bs[pop]["bs_within"]["synonymous"])
        err.append(1.96 * bs[pop]["bs_within"]["missense"])
        err.append(1.96 * bs[pop]["bs_between"]["synonymous"])
        err.append(1.96 * bs[pop]["bs_between"]["missense"])

yy = np.array(yy)

ax1.plot([0, 65], [0, 0], "k--", lw=lw)

for ii, v in enumerate(xx):
    if ii % 4 == 0:
        xx[ii] = v + 0.4
    elif ii % 4 == 1:
        xx[ii] = v + 0.2
    elif ii % 4 == 2:
        xx[ii] = v - 0.2
    elif ii % 4 == 3:
        xx[ii] = v - 0.4

ax1.plot(
    xx.reshape(15, 4)[:, 0],
    yy.reshape(15, 4)[:, 0],
    "o",
    ms=ms,
    color=colors[0],
    label="Synonymous, within domains",
)
ax1.plot(
    xx.reshape(15, 4)[:, 1],
    yy.reshape(15, 4)[:, 1],
    "s",
    ms=ms,
    color=colors[1],
    label="Missense, within",
)
ax1.plot(
    xx.reshape(15, 4)[:, 2],
    yy.reshape(15, 4)[:, 2],
    "o",
    ms=ms,
    markerfacecolor="white",
    color=colors[0],
    label="Synonymous, between",
)
ax1.plot(
    xx.reshape(15, 4)[:, 3],
    yy.reshape(15, 4)[:, 3],
    "s",
    ms=ms,
    markerfacecolor="white",
    color=colors[1],
    label="Missense, between",
)


ax1.errorbar(xx, yy, yerr=err, linestyle="None", linewidth=lw, color="gray")

ax1.set_ylabel(r"$\sigma_d^1$")
ax1.set_xlabel("Populations")

ax1.set_xticks(np.mean(xx.reshape(15, 4), axis=1))
ax1.set_xticklabels(
    populations["Africa"] + populations["Europe"] + populations["East Asia"]
)

# ax1.text(xx[0], 0.1, "synonymous", rotation=90, ha="center", va="center")
# ax1.text(xx[1], 0.12, "missense", rotation=90, ha="center", va="center")
# ax1.text(5, 0.005, "neutral expectation", ha="center", va="center")

ax1.legend(loc="upper left")

ax1.set_xlim([0, 65])
ax1.set_title("Signed LD within and between domains")

ax2 = plt.subplot2grid((5, 3), (3, 0), rowspan=2)

plot_inside_outside("YRI", ax2, ylabel=True, legend=True)

ax3 = plt.subplot2grid((5, 3), (3, 1), rowspan=2)

plot_inside_outside("IBS", ax3)

ax4 = plt.subplot2grid((5, 3), (3, 2), rowspan=2)

plot_inside_outside("KHV", ax4)

fig.tight_layout()
fig.text(0.02, 0.95, "A", fontsize=8, ha="center", va="center")
fig.text(0.05, 0.38, "B", fontsize=8, ha="center", va="center")
fig.text(0.38, 0.38, "C", fontsize=8, ha="center", va="center")
fig.text(0.70, 0.38, "D", fontsize=8, ha="center", va="center")

#plt.savefig(f"fig6.pdf")
plt.savefig(f"fig6.png", dpi=300)
