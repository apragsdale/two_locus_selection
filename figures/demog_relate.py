# Figure 4: Demography and trajectories.
# Two columns, left is simple size changes: expansion and bottleneck+recovery
# Right column has human-like domographic history inferred by Relate
# Rows: all with rho = 1
# 1: illustration of demography
# 2: (a) Additive selection w/ gamma=-1, (b) negative epistasis, gamma=-1, (c) negative epistasis, gamma=-10
# 3: (a) Simple dominance, gamma=-1, (b) simple dominance, gamma=-10
# 4: (a) Some gene-based selection model
import sys
import moments.TwoLocus
import matplotlib.pylab as plt
import numpy as np
import pickle
import copy

# set font sizes
import matplotlib

plt.rcParams["legend.title_fontsize"] = "xx-small"
matplotlib.rc("xtick", labelsize=6)
matplotlib.rc("ytick", labelsize=6)
matplotlib.rc("axes", labelsize=7)
matplotlib.rc("axes", titlesize=7)
matplotlib.rc("legend", fontsize=6)

from bokeh.palettes import Colorblind

n = 30  # the sample size


def piecewise_constant(
    nus,
    Ts,
    n,
    rho=0.0,
    gamma=0,
    h=None,
    dom="simple",
    eps=None,
    track_stats=True,
    spacing=3,
    conditions=[],
):
    if track_stats:
        spectra = {"t": [], "data": []}
    if gamma == 0:
        sel_params = None
        sel_params_general = None
    elif h is None and eps is None:
        sel_params = moments.TwoLocus.Util.additive_epistasis(gamma)
        sel_params_general = None
    elif h is not None:
        sel_params = None
        if dom == "simple":
            sel_params_general = moments.TwoLocus.Util.simple_dominance(gamma, h=h)
        elif dom == "gene":
            sel_params_general = moments.TwoLocus.Util.gene_based_dominance(gamma, h=h)
    elif eps is not None:
        sel_params = moments.TwoLocus.Util.additive_epistasis(gamma, epsilon=eps)
        sel_params_general = None
    F = moments.TwoLocus.Demographics.equilibrium(
        n, rho=rho, sel_params=sel_params, sel_params_general=sel_params_general
    )
    if track_stats:
        spectra["t"].append(0)
        spectra["data"].append(copy.deepcopy(F))

    print("o" * len(nus), flush=True)
    for nu, T in zip(nus, Ts):
        print(".", end="", flush=True)
        if not track_stats:
            F.integrate(
                nu,
                T,
                rho=rho,
                sel_params=sel_params,
                sel_params_general=sel_params_general,
            )
        else:
            T_sub = T / spacing
            for i in range(spacing):
                F.integrate(
                    nu,
                    T_sub,
                    rho=rho,
                    sel_params=sel_params,
                    sel_params_general=sel_params_general,
                )
                spectra["t"].append(spectra["t"][-1] + T_sub)
                spectra["data"].append(copy.deepcopy(F))
    print()
    return spectra


# load them:
def load_relate_curves(fname, pop0, pop1):
    t = [-1]
    ind0 = -1
    ind1 = -1
    for line in open(fname):
        if pop0 in line and pop1 in line:
            ind0 = line.split().index(pop0)
            ind1 = line.split().index(pop1)
        elif t[0] == -1 and ind0 != -1 and ind1 != -1:
            t = np.array([float(x) for x in line.split()])
        else:
            if line.startswith(f"{ind0} {ind0}"):
                coal0 = line.split()[2:]
                N0 = 1 / 2 / np.array([float(x) for x in coal0])
            elif line.startswith(f"{ind1} {ind1}"):
                coal1 = line.split()[2:]
                N1 = 1 / 2 / np.array([float(x) for x in coal1])
    return t, N0, N1


## Relate demography:
pop0 = "YRI"
pop1 = "CEU"
t, N0, N1 = load_relate_curves(
    "../../african-demography/coalescent_rate_curves/agr_refit_1_22.coal", pop0, pop1
)

# adjust to set recent size to reasonable value, set Ne for t > 1e6 years to

t = np.concatenate(([0], t[4:23]))
Ne = np.mean([N0[-1], N1[-1]])
N0 = N0[3:23]
N1 = N1[3:23]
Ne = np.mean([N0[-1], N1[-1]])
N0 = np.concatenate((N0, [Ne]))
N1 = np.concatenate((N1, [Ne]))

# convert to forward genetic units

Ts = ((t[1:] - t[:-1]) / 2 / Ne)[::-1]
nu0s = N0[::-1][2:] / Ne
nu1s = N1[::-1][2:] / Ne

## function to plot t vs N


def plot_relate_curve(ax, t, N, gen=29, line_style="--", lw=1.5, color="k", label=None):
    for ii, (l, r) in enumerate(zip(t[:-1], t[1:])):
        Ne = N[ii]
        if np.isinf(Ne):
            jj = 0
            Ne = N[jj]
            while np.isinf(Ne):
                jj += 1
                Ne = N[jj]
        if ii == 0:
            label = label
        else:
            label = None
        ax.plot(
            (l * gen, r * gen), (Ne, Ne), line_style, lw=lw, color=color, label=label
        )
        ax.plot(
            (r * gen, r * gen),
            (Ne, N[ii + 1]),
            line_style,
            lw=lw,
            color=color,
            label=None,
        )


### Compute data
data_fname = f"fig4_data_n_{n}.bp"
try:
    data = pickle.load(open(data_fname, "rb"))
except:
    data = {
        "Expansion": {"A": {}, "B": {}, "C": {}},
        "Bottleneck": {"A": {}, "B": {}, "C": {}},
        pop0: {"A": {}, "B": {}, "C": {}},
        pop1: {"A": {}, "B": {}, "C": {}},
    }

    ## second row (A): rho=0, rho=5 with gamma = -2
    print("Additive")
    # Expansion model
    print("Expansion")
    F, stats = piecewise_constant(nus_expand, Ts_expand, n, rho=0, gamma=-2, spacing=10)
    data["Expansion"]["A"][0] = stats
    F, stats = piecewise_constant(nus_expand, Ts_expand, n, rho=5, gamma=-2, spacing=10)
    data["Expansion"]["A"][5] = stats
    # Bottleneck model
    print("Bottleneck")
    F, stats = piecewise_constant(nus_bottle, Ts_bottle, n, rho=0, gamma=-2, spacing=10)
    data["Bottleneck"]["A"][0] = stats
    F, stats = piecewise_constant(nus_bottle, Ts_bottle, n, rho=5, gamma=-2, spacing=10)
    data["Bottleneck"]["A"][5] = stats
    # African model
    print(pop0)
    F, stats = piecewise_constant(nu0s, Ts, n, rho=0, gamma=-2, spacing=10)
    data[pop0]["A"][0] = stats
    F, stats = piecewise_constant(nu0s, Ts, n, rho=5, gamma=-2, spacing=10)
    data[pop0]["A"][5] = stats
    # European model
    print(pop1)
    F, stats = piecewise_constant(nu1s, Ts, n, rho=0, gamma=-2, spacing=10)
    data[pop1]["A"][0] = stats
    F, stats = piecewise_constant(nu1s, Ts, n, rho=5, gamma=-2, spacing=10)
    data[pop1]["A"][5] = stats

    ## Third row (B): rho=1, gamma = -2, eps=0.5 and -0.5
    print("B: epistasis")
    # Expansion model
    print("Expansion")
    F, stats = piecewise_constant(
        nus_expand, Ts_expand, n, rho=0, gamma=-2, eps=0.5, spacing=10
    )
    data["Expansion"]["B"][0.5] = stats
    F, stats = piecewise_constant(
        nus_expand, Ts_expand, n, rho=5, gamma=-2, eps=-0.5, spacing=10
    )
    data["Expansion"]["B"][-0.5] = stats
    # Bottleneck model
    print("Bottleneck")
    F, stats = piecewise_constant(
        nus_bottle, Ts_bottle, n, rho=0, gamma=-2, eps=0.5, spacing=10
    )
    data["Bottleneck"]["B"][0.5] = stats
    F, stats = piecewise_constant(
        nus_bottle, Ts_bottle, n, rho=5, gamma=-2, eps=-0.5, spacing=10
    )
    data["Bottleneck"]["B"][-0.5] = stats
    # African model
    print(pop0)
    F, stats = piecewise_constant(nu0s, Ts, n, rho=0, gamma=-2, eps=0.5, spacing=10)
    data[pop0]["B"][0.5] = stats
    F, stats = piecewise_constant(nu0s, Ts, n, rho=5, gamma=-2, eps=-0.5, spacing=10)
    data[pop0]["B"][-0.5] = stats
    # European model
    print(pop1)
    F, stats = piecewise_constant(nu1s, Ts, n, rho=0, gamma=-2, eps=0.5, spacing=10)
    data[pop1]["B"][0.5] = stats
    F, stats = piecewise_constant(nu1s, Ts, n, rho=5, gamma=-2, eps=-0.5, spacing=10)
    data[pop1]["B"][-0.5] = stats

    ## Fourth row (C): rho=1, gamma=-2, h=0.1, with (i) simple, and (ii) gene based dom.
    # Expansion model
    print("C: dominance")
    print("Expansion")
    F, stats = piecewise_constant(
        nus_expand, Ts_expand, n, rho=0, gamma=-2, h=0.1, dom="simple", spacing=10
    )
    data["Expansion"]["C"]["simple"] = stats
    F, stats = piecewise_constant(
        nus_expand, Ts_expand, n, rho=5, gamma=-2, h=0.1, dom="gene", spacing=10
    )
    data["Expansion"]["C"]["gene"] = stats
    # Bottleneck model
    print("Bottleneck")
    F, stats = piecewise_constant(
        nus_bottle, Ts_bottle, n, rho=0, gamma=-2, h=0.1, dom="simple", spacing=10
    )
    data["Bottleneck"]["C"]["simple"] = stats
    F, stats = piecewise_constant(
        nus_bottle, Ts_bottle, n, rho=5, gamma=-2, h=0.1, dom="gene", spacing=10
    )
    data["Bottleneck"]["C"]["gene"] = stats
    # African model
    print(pop0)
    F, stats = piecewise_constant(nu0s, Ts, n, rho=0, gamma=-2, h=0.1, dom="simple", spacing=10)
    data[pop0]["C"]["simple"] = stats
    F, stats = piecewise_constant(nu0s, Ts, n, rho=5, gamma=-2, h=0.1, dom="gene", spacing=10)
    data[pop0]["C"]["gene"] = stats
    # European model
    print(pop1)
    F, stats = piecewise_constant(nu1s, Ts, n, rho=0, gamma=-2, h=0.1, dom="simple", spacing=10)
    data[pop1]["C"]["simple"] = stats
    F, stats = piecewise_constant(nu1s, Ts, n, rho=5, gamma=-2, h=0.1, dom="gene", spacing=10)
    data[pop1]["C"]["gene"] = stats

fig = plt.figure(4, figsize=(6.5, 6))
fig.clf()

markers = ["x", "+", "1"]
colors = Colorblind[8]

ax1 = plt.subplot(4, 2, 1)

plot_relate_curve(
    ax1,
    t_expand,
    N_expand,
    line_style="-",
    lw=1,
    color=colors[0],
    label="Expansion",
    gen=1,
)
plot_relate_curve(
    ax1,
    t_bottle,
    N_bottle,
    line_style="--",
    lw=1,
    color=colors[1],
    label="Bottleneck",
    gen=1,
)

ax1.legend()
# ax1.set_xscale("log")
# ax1.set_yscale("log")
ax1.set_xlim(4000, 0)
ax1.set_ylim(0, 5e4)
ax1.set_ylabel(r"$N_e$")

ax2 = plt.subplot(4, 2, 2)

plot_relate_curve(ax2, t, N0, line_style="-", lw=1, color=colors[0], label=pop0)
plot_relate_curve(ax2, t, N1, line_style="--", lw=1, color=colors[1], label=pop1)

ax2.legend()
ax2.set_xscale("log")
ax2.set_yscale("log")
ax2.set_xlim(1e6, 1e3)
ax2.set_ylim(1e2, 1e7)
ax2.set_ylabel(r"$N_e$")

Ts_gens = {}
Ts_gens["Expansion"] = 4000 - 1e4 * np.array(data["Expansion"]["A"][0]["t"])
Ts_gens["Bottleneck"] = 4000 - 1e4 * np.array(data["Bottleneck"]["A"][0]["t"])
Ts_years = 1e6 - 2 * Ne * 29 * np.array(data[pop0]["A"][0]["t"])

ax3 = plt.subplot(4, 2, 3)

ax3.plot(Ts_gens["Expansion"], 0 * Ts_gens["Expansion"], "k:", lw=1, label=None)
for rho, c in zip([0, 5], [colors[2], colors[3]]):
    for pop, m in zip(["Expansion", "Bottleneck"], ["-", "--"]):
        ax3.plot(
            Ts_gens[pop],
            data[pop]["A"][rho]["sd1"]["all"],
            m,
            lw=1,
            color=c,
            label=rf"{pop[:3]}, $\rho={rho}$",
        )

ax3.legend(ncol=2)
ax3.set_xlim(ax1.get_xlim())
ax3.set_ylabel(f"$\sigma_d^1$")
ax3.set_title("Additive selection, no epistasis")


ax4 = plt.subplot(4, 2, 4)
ax4.plot(Ts_years, 0 * Ts_years, "k:", lw=1, label=None)
for rho, c in zip([0, 5], [colors[2], colors[3]]):
    for pop, m in zip([pop0, pop1], ["-", "--"]):
        ax4.plot(
            Ts_years,
            data[pop]["A"][rho]["sd1"]["all"],
            m,
            lw=1,
            color=c,
            label=rf"{pop[:3]}, $\rho={rho}$",
        )

ax4.set_xlim(ax2.get_xlim())
ax4.set_title("Additive selection, no epistasis")
ax4.set_xscale("log")
ax4.legend(ncol=2)

ax5 = plt.subplot(4, 2, 5)

ax5.plot(Ts_gens["Expansion"], 0 * Ts_gens["Expansion"], "k:", lw=1, label=None)
for eps, c in zip([-0.5, 0.5], [colors[4], colors[5]]):
    for pop, m in zip(["Expansion", "Bottleneck"], ["-", "--"]):
        ax5.plot(
            Ts_gens[pop],
            data[pop]["B"][eps]["sd1"]["all"],
            m,
            lw=1,
            color=c,
            label=rf"{pop[:3]}, $\epsilon={eps}$",
        )

ax5.legend(ncol=2)
ax5.set_xlim(ax1.get_xlim())
ax5.set_ylabel(f"$\sigma_d^1$")
ax5.set_title(r"Epistasis (with $\rho=1$)")


ax6 = plt.subplot(4, 2, 6)
ax6.plot(Ts_years, 0 * Ts_years, "k:", lw=1, label=None)
for eps, c in zip([-0.5, 0.5], [colors[4], colors[5]]):
    for pop, m in zip([pop0, pop1], ["-", "--"]):
        ax6.plot(
            Ts_years,
            data[pop]["B"][eps]["sd1"]["all"],
            m,
            lw=1,
            color=c,
            label=rf"{pop[:3]}, $\epsilon={eps}$",
        )

ax6.set_xlim(ax2.get_xlim())
ax6.set_title(r"Epistasis (with $\rho=1$)")
ax6.set_xscale("log")
ax6.legend(ncol=2)

ax7 = plt.subplot(4, 2, 7)

ax7.plot(Ts_gens["Expansion"], 0 * Ts_gens["Expansion"], "k:", lw=1, label=None)
for model, c in zip(["simple", "gene"], [colors[6], colors[7]]):
    if model == "simple":
        model_type = "within loci"
    elif model == "gene":
        model_type = "gene-based"
    for pop, m in zip(["Expansion", "Bottleneck"], ["-", "--"]):
        ax7.plot(
            Ts_gens[pop],
            data[pop]["C"][model]["sd1"]["all"],
            m,
            lw=1,
            color=c,
            label=rf"{pop[:3]}, {model_type}",
        )

ax7.legend(ncol=2)
ax7.set_xlim(ax1.get_xlim())
ax7.set_ylabel(f"$\sigma_d^1$")
ax7.set_title(r"Dominance (with $\rho=1$)")
ax7.set_xlabel("Generations ago")

ax8 = plt.subplot(4, 2, 8)
ax8.plot(Ts_years, 0 * Ts_years, "k:", lw=1, label=None)
for model, c in zip(["simple", "gene"], [colors[6], colors[7]]):
    if model == "simple":
        model_type = "within loci"
    elif model == "gene":
        model_type = "gene-based"
    for pop, m in zip([pop0, pop1], ["-", "--"]):
        ax8.plot(
            Ts_years,
            data[pop]["C"][model]["sd1"]["all"],
            m,
            lw=1,
            color=c,
            label=rf"{pop[:3]}, {model_type}",
        )

ax8.set_xlim(ax2.get_xlim())
ax8.set_title(r"Dominance (with $\rho=1$)")
ax8.set_xscale("log")
ax8.set_xlabel("Time ago (years)")
ax8.legend(ncol=2)

fig.tight_layout()
fig.text(0.02, 0.98, "A", fontsize=8, ha="center", va="center")
fig.text(0.52, 0.98, "B", fontsize=8, ha="center", va="center")
fig.text(0.02, 0.73, "C", fontsize=8, ha="center", va="center")
fig.text(0.52, 0.73, "D", fontsize=8, ha="center", va="center")
fig.text(0.02, 0.48, "E", fontsize=8, ha="center", va="center")
fig.text(0.52, 0.48, "F", fontsize=8, ha="center", va="center")
fig.text(0.02, 0.23, "G", fontsize=8, ha="center", va="center")
fig.text(0.52, 0.23, "H", fontsize=8, ha="center", va="center")

plt.savefig(f"fig4_n0_{n}_toy_relate.pdf")
