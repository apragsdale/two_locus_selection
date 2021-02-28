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

fig = plt.figure(5, figsize=(6.5, 4.5))
fig.clf()

markers = ["x", "+", "1"]
colors = Colorblind[8]

# plot all stats within genes
ax1 = plt.subplot2grid((7, 4), (0, 0), rowspan=2, colspan=4)
ax1b = plt.subplot2grid((7, 4), (2, 0), rowspan=2, colspan=4)

xx = np.concatenate(
    (np.linspace(1, 15, 15), np.linspace(18, 32, 15), np.linspace(35, 49, 15))
)
yy = []
err = []

for continent in ["Africa", "Europe", "East Asia"]:
    for pop in populations[continent]:
        yy.append(
            data[pop]["ld_within"]["synonymous"][3]
            / data[pop]["ld_within"]["synonymous"][2]
        )
        yy.append(
            data[pop]["ld_within"]["missense"][3]
            / data[pop]["ld_within"]["missense"][2]
        )
        yy.append(
            data[pop]["ld_within"]["loss_of_function"][3]
            / data[pop]["ld_within"]["loss_of_function"][2]
        )
        err.append(1.96 * bs[pop]["bs_within"]["synonymous"])
        err.append(1.96 * bs[pop]["bs_within"]["missense"])
        err.append(1.96 * bs[pop]["bs_within"]["loss_of_function"])

yy = np.array(yy)
err = np.array(err)

ax1.plot([0, 50], [0, 0], "k--", lw=lw)
ax1b.plot([0, 50], [0, 0], "k--", lw=lw)

for ii, v in enumerate(xx):
    if ii % 3 == 0:
        xx[ii] = v + 0.25
    elif ii % 3 == 2:
        xx[ii] = v - 0.25

ax1.plot(
    xx.reshape(15, 3)[:, 0],
    yy.reshape(15, 3)[:, 0],
    "o",
    ms=ms,
    color=colors[0],
    label="Synonymous",
)
ax1.plot(
    xx.reshape(15, 3)[:, 1],
    yy.reshape(15, 3)[:, 1],
    "s",
    ms=ms,
    color=colors[1],
    label="Missense",
)
ax1b.plot(
    xx.reshape(15, 3)[:, 2],
    yy.reshape(15, 3)[:, 2],
    "X",
    ms=ms,
    color=colors[3],
    label="Loss of function",
)

ax1.errorbar(
    xx.reshape(15, 3)[:, 0],
    yy.reshape(15, 3)[:, 0],
    yerr=err.reshape(15, 3)[:, 0],
    linestyle="None",
    linewidth=lw,
    color="gray",
)
ax1.errorbar(
    xx.reshape(15, 3)[:, 1],
    yy.reshape(15, 3)[:, 1],
    yerr=err.reshape(15, 3)[:, 1],
    linestyle="None",
    linewidth=lw,
    color="gray",
)
ax1b.errorbar(
    xx.reshape(15, 3)[:, 2],
    yy.reshape(15, 3)[:, 2],
    yerr=err.reshape(15, 3)[:, 2],
    linestyle="None",
    linewidth=lw,
    color="gray",
)

ax1.set_ylabel(r"$\sigma_d^1$")
ax1b.set_ylabel(r"$\sigma_d^1$")
ax1b.set_xlabel("Populations")

ax1.set_xticks(
    [2, 5, 8, 11, 14, 19, 22, 25, 28, 31, 36, 39, 42, 45, 48,]
)
ax1.set_xticklabels([])
ax1b.set_xticks(
    [2, 5, 8, 11, 14, 19, 22, 25, 28, 31, 36, 39, 42, 45, 48,]
)
ax1b.set_xticklabels(
    populations["Africa"] + populations["Europe"] + populations["East Asia"]
)

# ax1.text(xx[0], 0.1, "synonymous", rotation=90, ha="center", va="center")
# ax1.text(xx[1], 0.12, "missense", rotation=90, ha="center", va="center")
# ax1.text(5, 0.005, "neutral expectation", ha="center", va="center")

ax1.legend(loc="upper left")
ax1b.legend(loc="upper left")

ax1.set_xlim([0, 50])
ax1b.set_xlim(ax1.get_xlim())
ax1.set_title("Signed LD within genes")
# ax1.set_ylim(-0.9, 0.2)

## plots between genes


def plot_between(pop, data, bs, ax, legend=False, ylabel=False):
    for annot, marker, color in zip(
        ["synonymous", "missense"], ["o", "s", "X"], colors[:3]
    ):
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
        jitter = 0.1 - 0.2 * (annot == "synonymous")
        annot = annot[0].upper() + annot[1:]
        ax.plot(
            xx + jitter, to_plot, marker + "-", ms=ms, color=color, lw=lw, label=annot
        )
        ax.errorbar(
            xx + jitter, to_plot, yerr=err, linestyle="None", linewidth=lw, color="gray"
        )
    ax.set_xticks(xx)
    ax.set_xticklabels(labels, rotation=45, ha="right")
    ax.set_xlabel("Distance between bp")
    if legend:
        ax.legend(loc="upper right")
    ax.plot(xx, np.zeros(len(xx)), "k--", lw=lw, label=None)
    if ylabel:
        ax.set_ylabel("$\sigma_d^1$")
    ax.set_title(pop)


ax2 = plt.subplot2grid((7, 4), (5, 0), rowspan=2)
pop = "ESN"
plot_between(pop, data, bs, ax2, legend=True, ylabel=True)

ax3 = plt.subplot2grid((7, 4), (5, 1), rowspan=2)
pop = "MSL"
plot_between(pop, data, bs, ax3)

ax4 = plt.subplot2grid((7, 4), (5, 2), rowspan=2)
pop = "GBR"
plot_between(pop, data, bs, ax4)

ax5 = plt.subplot2grid((7, 4), (5, 3), rowspan=2)
pop = "CHB"
plot_between(pop, data, bs, ax5)

fig.tight_layout()
fig.subplots_adjust(hspace=0.3)
plt.savefig(f"fig5.pdf")
