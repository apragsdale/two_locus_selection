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
data_dom = {}
data_btw = {}
SE = {}
SE_dom = {}
SE_btw = {}

populations = {
    "Africa": ["ESN", "GWD", "LWK", "MSL", "YRI"],
    "Europe": ["CEU", "FIN", "GBR", "IBS", "TSI"],
    "East Asia": ["CDX", "CHB", "CHS", "JPT", "KHV"],
}

for continent, pops in populations.items():
    for pop in pops:
        ld_pop = pickle.load(
            open(f"../analysis/parsed_data_v2/{pop}.unphased.all.bp", "rb")
        )
        data[pop] = ld_pop
        se_pop = pickle.load(
            open(
                f"../analysis/parsed_data_v2/{pop}.unphased.all.SEs.bp",
                "rb",
            )
        )
        SE[pop] = se_pop

        ld_pop = pickle.load(
            open(f"../analysis/parsed_data_v2/{pop}.unphased.domains.bp", "rb")
        )
        data_dom[pop] = ld_pop
        se_pop = pickle.load(
            open(
                f"../analysis/parsed_data_v2/{pop}.unphased.domains.SEs.bp",
                "rb",
            )
        )
        SE_dom[pop] = se_pop

        ld_pop = pickle.load(
            open(f"../analysis/parsed_data_v2/{pop}.unphased.between_domains.bp", "rb")
        )
        data_btw[pop] = ld_pop
        se_pop = pickle.load(
            open(
                f"../analysis/parsed_data_v2/{pop}.unphased.between_domains.SEs.bp",
                "rb",
            )
        )
        SE_btw[pop] = se_pop

ms = 2
lw = 0.5
markers = ["o", "s"]
colors = Colorblind[8]
labels = [
    "Syn. in domains",
    "Mis. in domains",
    "Syn. outside",
    "Mis. outside",
]


def plot_domain_data(pops, ax, legend=False):
    categ = "all"
    # within syn, within non, match syn, match non
    xx = np.array([1, 2, 3, 4])
    for i in range(len(pops) - 1):
        xx = np.concatenate((xx, xx[-1] + 1 + np.array([1, 2, 3, 4])))

    to_plot = [
        [
            data_dom[pop]["within"][categ]["synonymous"]["sd1"],
            data_dom[pop]["within"][categ]["missense"]["sd1"],
            data_btw[pop]["between"][categ]["synonymous"]["sd1"],
            data_btw[pop]["between"][categ]["missense"]["sd1"],
        ]
        for pop in pops
    ]
    err_plot = [
        [
            SE_dom[pop]["within"][categ]["synonymous"]["sd1"],
            SE_dom[pop]["within"][categ]["missense"]["sd1"],
            SE_btw[pop]["between"][categ]["synonymous"]["sd1"],
            SE_btw[pop]["between"][categ]["missense"]["sd1"],
        ]
        for pop in pops
    ]
    to_plot = np.array(to_plot).flatten()
    err_plot = np.array(err_plot).flatten()

    ax.plot(
        (xx[0] - 1, xx[-1] + 1),
        (0, 0),
        "k--",
        lw=lw,
    )

    for ii, (x, p, e) in enumerate(zip(xx, to_plot, err_plot)):
        if (ii / 4) % 1 >= 0.5:
            mfc = "white"
        else:
            mfc = None
        if ii < 4:
            label = labels[ii]
        else:
            label = None
        ax.plot(
            x,
            p,
            markers[ii % 2],
            ms=ms,
            color=colors[ii % 2],
            label=label,
            mfc=mfc,
        )
        ax.errorbar(
            x,
            p,
            yerr=e,
            linestyle="None",
            linewidth=lw,
            color="gray",
        )

    ax.set_xlim(xx[0] - 1, xx[-1] + 1)
    ax.set_xticks(np.mean(xx.reshape(len(pops), 4), axis=1))
    ax.set_xticklabels(pops)
    if legend:
        ax.legend(fontsize=5, frameon=False, loc="upper left")


## unplotted populations, domain info
fig = plt.figure(849275, figsize=(8, 5))
fig.clf()
i = 0
for cont in ["Africa", "Europe", "East Asia"]:
    for pop in populations[cont]:
        ax = plt.subplot(3, 5, i + 1)
        ax.set_title(pop)
        plot_domain_data([pop], ax)
        if i % 5 == 0:
            ax.set_ylabel(r"$\sigma_d^1$")
        ax.set_xticks([1, 2, 3, 4])
        if i < 10:
            ax.set_xticklabels([])
        else:
            ax.set_xticklabels(["Syn.\nw/in", "Mis.\nw/in", "Syn.\nbtw.", "Mis.\nbtw."])
        i += 1


plt.tight_layout()
plt.savefig("data_within_between.pdf")
