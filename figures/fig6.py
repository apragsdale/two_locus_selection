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

fig = plt.figure(5, figsize=(6.5, 3))
fig.clf()

markers = ["x", "+", "1"]
colors = Colorblind[8]

# plot all stats within genes
ax1 = plt.subplot(1, 1, 1)

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

fig.tight_layout()

plt.savefig(f"fig6.pdf")
