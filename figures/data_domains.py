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

SE = {}

populations = {
    "Africa": ["ESN", "GWD", "LWK", "MSL", "YRI"],
    "Europe": ["CEU", "FIN", "GBR", "IBS", "TSI"],
    "East Asia": ["CDX", "CHB", "CHS", "JPT", "KHV"],
}

for continent, pops in populations.items():
    for pop in pops:
        ld_pop = pickle.load(
            open(f"../analysis/parsed_data_v2/{pop}.unphased.domains.bp", "rb")
        )
        data[pop] = ld_pop
        se_pop = pickle.load(
            open(
                f"../analysis/parsed_data_v2/{pop}.unphased.domains.SEs.bp",
                "rb",
            )
        )
        SE[pop] = se_pop

ms = 3
lw = 0.5
markers = ["o", "s"]
colors = Colorblind[8]
labels = [
    "Synonymous, w/in domains",
    "Missense, w/in domains",
    "Syn., outside domains",
    "Mis., outside domains",
]


def plot_data(pop, categ, ax):
    # within syn, within non, match syn, match non
    to_plot = [
        data[pop]["within"][categ]["synonymous"]["sd1"],
        data[pop]["within"][categ]["missense"]["sd1"],
        data[pop]["outside"]["matching"][categ]["synonymous"]["sd1"],
        data[pop]["outside"]["matching"][categ]["missense"]["sd1"],
    ]
    err_plot = [
        SE[pop]["within"][categ]["synonymous"]["sd1"],
        SE[pop]["within"][categ]["missense"]["sd1"],
        SE[pop]["outside"]["matching"][categ]["synonymous"]["sd1"],
        SE[pop]["outside"]["matching"][categ]["missense"]["sd1"],
    ]

    for ii, (p, e) in enumerate(zip(to_plot, err_plot)):
        if ii >= 2:
            mfc = "white"
        else:
            mfc = None
        ax.plot(
            ii,
            p,
            ms=ms,
            color=colors[ii % 2],
            marker=markers[ii % 2],
            label=labels[ii],
            mfc=mfc,
        )
        ax.errorbar(
            ii,
            p,
            yerr=e,
            linestyle="None",
            linewidth=lw,
            color="gray",
        )
    ax.hlines(
        0, ax.get_xlim()[0], ax.get_xlim()[1], linestyle="--", colors=["black"], lw=1
    )


fig = plt.figure(9876, figsize=(6.5, 4.5))
plt.clf()
pop_afr = "ESN"
pop_eur = "CEU"
pop_eas = "CDX"

ax1 = plt.subplot(3, 4, 1)
categ = "all"
plot_data(pop_afr, categ, ax1)
ax1.set_ylabel(f"{pop_afr}\n" + r"$\sigma_d^1$")
ax1.set_title("All pairs")

ax2 = plt.subplot(3, 4, 2)
categ = "leq2"
plot_data(pop_afr, categ, ax2)
ax2.set_title(r"$n_A, n_B \leq 2$")

ax3 = plt.subplot(3, 4, 3)
categ = "3to8"
plot_data(pop_afr, categ, ax3)
ax3.set_title(r"$3 \leq n_A, n_B \leq 8$")

ax4 = plt.subplot(3, 4, 4)
categ = "geq9"
plot_data(pop_afr, categ, ax4)
ax4.set_title(r"$n_A, n_B \geq 9$")

ax5 = plt.subplot(3, 4, 5)
categ = "all"
plot_data(pop_eur, categ, ax5)
ax5.set_ylabel(f"{pop_eur}\n" + r"$\sigma_d^1$")

ax6 = plt.subplot(3, 4, 6)
categ = "leq2"
plot_data(pop_eur, categ, ax6)

ax7 = plt.subplot(3, 4, 7)
categ = "3to8"
plot_data(pop_eur, categ, ax7)

ax8 = plt.subplot(3, 4, 8)
categ = "geq9"
plot_data(pop_eur, categ, ax8)

ax9 = plt.subplot(3, 4, 9)
categ = "all"
plot_data(pop_eas, categ, ax9)
ax9.set_ylabel(f"{pop_eas}\n" + r"$\sigma_d^1$")

ax10 = plt.subplot(3, 4, 10)
categ = "leq2"
plot_data(pop_eas, categ, ax10)

ax11 = plt.subplot(3, 4, 11)
categ = "3to8"
plot_data(pop_eas, categ, ax11)

ax12 = plt.subplot(3, 4, 12)
categ = "geq9"
plot_data(pop_eas, categ, ax12)


plt.tight_layout()
plt.savefig("data_domains.pdf")

fig2 = plt.figure(98765, figsize=(6.5, 7))
fig2.clf()
i = 0
for pop in populations["Africa"]:
    #if pop != pop_afr:
        for categ in ["all", "leq2", "3to8", "geq9"]:
            ax = plt.subplot(5, 4, i + 1)
            if i == 0:
                ax.set_title("All pairs")
            elif i == 1:
                ax.set_title(r"$n_A, n_B \leq 2$")
            elif i == 2:
                ax.set_title(r"$3 \leq n_A, n_B \leq 8$")
            elif i == 3:
                ax.set_title(r"$n_A, n_B \geq 9$")
            plot_data(pop, categ, ax)
            if i % 4 == 0:
                ax.set_ylabel(f"{pop}\n" + r"$\sigma_d^1$")
            ax.set_xticks([0, 1, 2, 3])
            if i < 12:
                ax.set_xticklabels([])
            else:
                ax.set_xticklabels(["Syn.\nw/in", "Mis.\nw/in", "Syn.\nout.", "Mis.\nout."])
            i += 1


plt.tight_layout()
plt.savefig("data_domains_afr.pdf")

fig3 = plt.figure(98764, figsize=(6.5, 7))
fig3.clf()
i = 0
for pop in populations["Europe"]:
    #if pop != pop_eur:
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
            plot_data(pop, categ, ax)
            if i % 4 == 0:
                ax.set_ylabel(f"{pop}\n" + r"$\sigma_d^1$")
            ax.set_xticks([0, 1, 2, 3])
            if i < 12:
                ax.set_xticklabels([])
            else:
                ax.set_xticklabels(["Syn.\nw/in", "Mis.\nw/in", "Syn.\nout.", "Mis.\nout."])
            i += 1

plt.tight_layout()
plt.savefig("data_domains_eur.pdf")

fig4 = plt.figure(98763, figsize=(6.5, 7))
fig4.clf()
i = 0
for pop in populations["East Asia"]:
    #if pop != pop_eas:
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
            plot_data(pop, categ, ax)
            if i % 4 == 0:
                ax.set_ylabel(f"{pop}\n" + r"$\sigma_d^1$")
            ax.set_xticks([0, 1, 2, 3])
            if i < 12:
                ax.set_xticklabels([])
            else:
                ax.set_xticklabels(["Syn.\nw/in", "Mis.\nw/in", "Syn.\nout.", "Mis.\nout."])
            i += 1

plt.tight_layout()
plt.savefig("data_domains_eas.pdf")

