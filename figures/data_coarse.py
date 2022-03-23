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
SE = {}
SE_dom = {}

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


ms = 3
lw = 0.5
markers = ["o", "s"]
colors = Colorblind[8]
labels = [
    "Syn. in domains",
    "Mis. in domains",
    "Syn. outside",
    "Mis. outside",
]


def plot_decay(pop, data, bs, ax, legend=False, ylabel=False):
    bins = data[pop]["bin_edges"][:-1]
    xx = np.mean(bins, axis=1)
    ax.plot(xx, np.zeros(len(xx)), "k--", lw=lw, label=None)
    for annot, marker, color in zip(["synonymous", "missense"], ["o", "s"], colors[:2]):
        data_pop = data[pop]["bins"]["all_sites"][annot]["sd1"][:-1]
        err_pop = [SE[pop]["bins"]["all_sites"][annot][b]["sd1"] for b in bins]
        annot = annot[0].upper() + annot[1:]
        ax.plot(xx, data_pop, marker + "-", ms=ms, color=color, lw=lw, label=annot)
        ax.errorbar(
            xx, data_pop, yerr=err_pop, linestyle="None", linewidth=lw, color="gray"
        )
    ax.set_xscale("log")
    ax.set_xlabel("Distance (bp)")
    if legend:
        ax.legend(loc="upper right", fontsize=5)
    if ylabel:
        ax.set_ylabel("$\sigma_d^1$")
    ax.set_title(pop)


def plot_domain_data(pops, categ, ax, legend=False):
    # within syn, within non, match syn, match non
    xx = np.array([1, 2, 3, 4])
    for i in range(len(pops) - 1):
        xx = np.concatenate((xx, xx[-1] + 1 + np.array([1, 2, 3, 4])))

    to_plot = [
        [
            data_dom[pop]["within"][categ]["synonymous"]["sd1"],
            data_dom[pop]["within"][categ]["missense"]["sd1"],
            data_dom[pop]["outside"]["matching"][categ]["synonymous"]["sd1"],
            data_dom[pop]["outside"]["matching"][categ]["missense"]["sd1"],
        ]
        for pop in pops
    ]
    err_plot = [
        [
            SE_dom[pop]["within"][categ]["synonymous"]["sd1"],
            SE_dom[pop]["within"][categ]["missense"]["sd1"],
            SE_dom[pop]["outside"]["matching"][categ]["synonymous"]["sd1"],
            SE_dom[pop]["outside"]["matching"][categ]["missense"]["sd1"],
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
        ax.legend(fontsize=5, frameon=False)


fig = plt.figure(5, figsize=(6.5, 4.5))
fig.clf()

# plot all stats within genes
dims = (8, 4)
ax1 = plt.subplot2grid(dims, (0, 0), rowspan=2, colspan=4)
ax1b = plt.subplot2grid(dims, (2, 0), rowspan=2, colspan=4)

xx = np.concatenate(
    (np.linspace(1, 15, 15), np.linspace(18, 32, 15), np.linspace(35, 49, 15))
)
yy = []
err = []

for continent in ["Africa", "Europe", "East Asia"]:
    for pop in populations[continent]:
        yy.append(data[pop]["all_sites"]["all"]["synonymous"]["sd1"])
        yy.append(data[pop]["all_sites"]["all"]["missense"]["sd1"])
        yy.append(data[pop]["all_sites"]["all"]["loss_of_function"]["sd1"])
        err.append(SE[pop]["all_sites"]["all"]["synonymous"]["sd1"])
        err.append(SE[pop]["all_sites"]["all"]["missense"]["sd1"])
        err.append(SE[pop]["all_sites"]["all"]["loss_of_function"]["sd1"])

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
    "D",
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
    [
        2,
        5,
        8,
        11,
        14,
        19,
        22,
        25,
        28,
        31,
        36,
        39,
        42,
        45,
        48,
    ]
)
ax1.set_xticklabels([])
ax1b.set_xticks(
    [
        2,
        5,
        8,
        11,
        14,
        19,
        22,
        25,
        28,
        31,
        36,
        39,
        42,
        45,
        48,
    ]
)
ax1b.set_xticklabels(
    populations["Africa"] + populations["Europe"] + populations["East Asia"]
)

# ax1.text(xx[0], 0.1, "synonymous", rotation=90, ha="center", va="center")
# ax1.text(xx[1], 0.12, "missense", rotation=90, ha="center", va="center")
# ax1.text(5, 0.005, "neutral expectation", ha="center", va="center")

ax1.legend(loc="upper left", fontsize=5)
ax1b.legend(loc="upper left", fontsize=5)

ax1.set_xlim([0, 50])
ax1b.set_xlim(ax1.get_xlim())
ax1.set_title("Signed LD within genes")
# ax1.set_ylim(-0.9, 0.2)

## plots of domain data
pop_afr = "YRI"
pop_eur = "CEU"
pop_eas = "JPT"
pops = [pop_afr, pop_eur, pop_eas]
ax2 = plt.subplot2grid(dims, (5, 0), rowspan=3)
plot_domain_data(pops, "all", ax2, legend=True)
ax2.set_ylabel("$\sigma_d^1$")
ax2.set_xlabel("Populations")
ax2.set_title("All pairs")

ax3 = plt.subplot2grid(dims, (5, 1), rowspan=3)
plot_domain_data(pops, "leq2", ax3)
ax3.set_xlabel("Populations")
ax3.set_title("$n_A, n_B \leq 2$")

ax4 = plt.subplot2grid(dims, (5, 2), rowspan=3)
plot_domain_data(pops, "3to8", ax4)
ax4.set_xlabel("Populations")
ax4.set_title("$3 \leq n_A, n_B \leq 8$")

ax5 = plt.subplot2grid(dims, (5, 3), rowspan=3)
plot_domain_data(pops, "geq9", ax5)
ax5.set_xlabel("Populations")
ax5.set_title("$n_A, n_B \geq 9$")

fig.tight_layout()
fig.subplots_adjust(hspace=0.5)
# fig.text(0.02, 0.95, "A", fontsize=8, ha="center", va="center")
# fig.text(0.02, 0.67, "B", fontsize=8, ha="center", va="center")
# fig.text(0.06, 0.33, "C", fontsize=8, ha="center", va="center")
# fig.text(0.30, 0.33, "D", fontsize=8, ha="center", va="center")
# fig.text(0.54, 0.33, "E", fontsize=8, ha="center", va="center")
# fig.text(0.77, 0.33, "F", fontsize=8, ha="center", va="center")

plt.savefig(f"data_compact.pdf")



## unplotted populations, domain info
fig2 = plt.figure(98765, figsize=(6.5, 6))
fig2.clf()
i = 0
for pop in populations["Africa"]:
    if pop != pop_afr:
        for categ in ["all", "leq2", "3to8", "geq9"]:
            ax = plt.subplot(4, 4, i + 1)
            if i == 0:
                ax.set_title("All pairs")
            elif i == 1:
                ax.set_title(r"$n_A, n_B \leq 2$")
            elif i == 2:
                ax.set_title(r"$3 \leq n_A, n_B \leq 8$")
            elif i == 3:
                ax.set_title(r"$n_A, n_B \geq 9$")
            plot_domain_data([pop], categ, ax)
            if i % 4 == 0:
                ax.set_ylabel(f"{pop}\n" + r"$\sigma_d^1$")
            ax.set_xticks([1, 2, 3, 4])
            if i < 12:
                ax.set_xticklabels([])
            else:
                ax.set_xticklabels(["Syn.\nw/in", "Mis.\nw/in", "Syn.\nout.", "Mis.\nout."])
            i += 1


plt.tight_layout()
plt.savefig("data_domains_afr.pdf")

fig3 = plt.figure(98764, figsize=(6.5, 6))
fig3.clf()
i = 0
for pop in populations["Europe"]:
    if pop != pop_eur:
        for categ in ["all", "leq2", "3to8", "geq9"]:
            ax = plt.subplot(4, 4, i + 1)
            if i == 0:
                ax.set_title("All pairs")
            elif i == 1:
                ax.set_title(r"$n_A, n_B \leq 2$")
            elif i == 2:
                ax.set_title(r"$3 \leq n_A, n_B \leq 8$")
            elif i == 3:
                ax.set_title(r"$n_A, n_B \geq 9$")
            plot_domain_data([pop], categ, ax)
            if i % 4 == 0:
                ax.set_ylabel(f"{pop}\n" + r"$\sigma_d^1$")
            ax.set_xticks([1, 2, 3, 4])
            if i < 12:
                ax.set_xticklabels([])
            else:
                ax.set_xticklabels(["Syn.\nw/in", "Mis.\nw/in", "Syn.\nout.", "Mis.\nout."])
            i += 1

plt.tight_layout()
plt.savefig("data_domains_eur.pdf")

fig4 = plt.figure(98763, figsize=(6.5, 6))
fig4.clf()
i = 0
for pop in populations["East Asia"]:
    if pop != pop_eas:
        for categ in ["all", "leq2", "3to8", "geq9"]:
            ax = plt.subplot(4, 4, i + 1)
            if i == 0:
                ax.set_title("All pairs")
            elif i == 1:
                ax.set_title(r"$n_A, n_B \leq 2$")
            elif i == 2:
                ax.set_title(r"$3 \leq n_A, n_B \leq 8$")
            elif i == 3:
                ax.set_title(r"$n_A, n_B \geq 9$")
            plot_domain_data([pop], categ, ax)
            if i % 4 == 0:
                ax.set_ylabel(f"{pop}\n" + r"$\sigma_d^1$")
            ax.set_xticks([1, 2, 3, 4])
            if i < 12:
                ax.set_xticklabels([])
            else:
                ax.set_xticklabels(["Syn.\nw/in", "Mis.\nw/in", "Syn.\nout.", "Mis.\nout."])
            i += 1

plt.tight_layout()
plt.savefig("data_domains_eas.pdf")

