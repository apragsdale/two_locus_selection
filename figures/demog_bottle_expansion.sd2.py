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
import gzip
import time

# set font sizes
import matplotlib

plt.rcParams["legend.title_fontsize"] = "xx-small"
matplotlib.rc("xtick", labelsize=6)
matplotlib.rc("ytick", labelsize=6)
matplotlib.rc("axes", labelsize=7)
matplotlib.rc("axes", titlesize=7)
matplotlib.rc("legend", fontsize=6)

from bokeh.palettes import Colorblind, Category10
from bokeh.palettes import YlGnBu, YlOrBr
c1 = YlGnBu[8]
c2 = YlOrBr[8]


def piecewise_constant(
    nus,
    Ts,
    n,
    rho=0.0,
    gamma=0,
    h=0.5,
    dom="simple",
    eps=0.0,
    spacing=3,
    name=None,
):
    if name is None:
        raise ValueError("need a name: bottle, expand, CEU, or YRI")

    fname = f"data/time_series/{name}-n_{n}-rho_{rho}-gamma_{gamma}-h_{h}-eps_{eps}-dom_{dom}-spacing_{spacing}.bp.gz"
    try:
        spectra = pickle.load(gzip.open(fname, "rb"))
        print("loaded from cached data")
    except IOError:
        spectra = {"t": [], "data": []}
        # set up selection parameters
        if gamma == 0:
            sel_params = None
            sel_params_general = None
        elif h != 0.5 and eps != 0.0:
            raise ValueError("can only have h or eps non-additive")
        elif h != 0.5:
            assert eps == 0.0
            sel_params = None
            if dom == "simple":
                sel_params_general = moments.TwoLocus.Util.simple_dominance(gamma, h=h)
            elif dom == "gene":
                sel_params_general = moments.TwoLocus.Util.gene_based_dominance(
                    gamma, h=h
                )
        else:
            sel_params = moments.TwoLocus.Util.additive_epistasis(gamma, epsilon=eps)
            sel_params_general = None
        # simulate data
        print(f"n={n}, rho={rho}, gamma={gamma}, h={h}, eps={eps}")
        print("o" * (len(nus) + 1), flush=True)
        F = moments.TwoLocus.Demographics.equilibrium(
            n, rho=rho, sel_params=sel_params, sel_params_general=sel_params_general
        )
        print(".", end="", flush=True)
        spectra["t"].append(0)
        spectra["data"].append(copy.deepcopy(F))

        for ii, (nu, T) in enumerate(zip(nus, Ts)):
            T_sub = T / spacing
            for i in range(spacing):
                if not (nu == 1 and np.all(np.array(nus[:ii]) == 1)):
                    F.integrate(
                        nu,
                        T_sub,
                        rho=rho,
                        sel_params=sel_params,
                        sel_params_general=sel_params_general,
                    )
                spectra["t"].append(spectra["t"][-1] + T_sub)
                spectra["data"].append(copy.deepcopy(F))
            print(".", end="", flush=True)
        print()
        pickle.dump(spectra, gzip.open(fname, "wb+"))
    return spectra


def ts_gens(t):
    return 1e4 * np.array(t)


## function to plot t vs N
def plot_relate_curve(ax, t, N, gen=29, line_style="--", lw=0.5, color="k", label=None):
    for ii, (l, r) in enumerate(zip(t[:-1], t[1:])):
        Ne = N[ii]
        # if np.isinf(Ne):
        #    jj = 0
        #    Ne = N[jj]
        #    while np.isinf(Ne):
        #        jj += 1
        #        Ne = N[jj]
        if ii == 0:
            label = label
        else:
            label = None
        ax.plot(
            (l * gen, r * gen), (Ne, Ne), line_style, lw=lw, color=color, label=label
        )
        if ii < len(t) - 2:
            ax.plot(
                (r * gen, r * gen),
                (Ne, N[ii + 1]),
                line_style,
                lw=lw,
                color=color,
                label=None,
            )


def get_sd2(
    nus, Ts, n, n_large=None, rho=0.0, gamma=0, eps=0.0, h=0.5, spacing=10, name=None
):
    if name is None:
        raise ValueError
    if n_large is not None:
        try:
            fname = f"data/time_series/{name}-n_{n}-rho_{rho}-gamma_{gamma}-h_{h}-eps_{eps}-dom_simple-spacing_{spacing}-large_{n_large}.bp.gz"
            data = pickle.load(gzip.open(fname))
        except IOError:
            print("could not find", fname)
            try:
                fname_large = f"data/time_series/{name}-n_{n_large}-rho_{rho}-gamma_{gamma}-h_{h}-eps_{eps}-dom_simple-spacing_{spacing}.bp.gz"
                data_large = pickle.load(gzip.open(fname_large, "rb"))
            except IOError:
                print("could not find", fname_large)
                print("calculating large")
                data_large = piecewise_constant(
                    nus,
                    Ts,
                    n_large,
                    gamma=gamma,
                    eps=eps,
                    h=h,
                    rho=rho,
                    spacing=spacing,
                    name=name,
                )
            print("Projecting")
            data = {"t": data_large["t"], "data": []}
            for F in data_large["data"]:
                print(" in ", end="", flush=True)
                time1 = time.time()
                data["data"].append(F.project(n))
                print(str(int(time.time() - time1)) + " seconds", flush=True)
            pickle.dump(data, gzip.open(fname, "wb+"))
    else:
        data = piecewise_constant(
            nus, Ts, n, gamma=gamma, eps=eps, h=h, rho=rho, spacing=spacing, name=name
        )
    v = []
    for F in data["data"]:
        v.append(F.D2() / F.pi2())
    ts = ts_gens(data["t"])
    return ts, v


def get_sd2_cond(
    rho=0.0,
    gamma=0,
    eps=0.0,
    h=0.5,
    spacing=10,
    name=None,
    nAmin=None,
    nBmin=None,
    nAmax=None,
    nBmax=None,
):
    n = 50
    data = pickle.load(
        gzip.open(
            f"data/time_series/{name}-n_{n}-rho_{rho}-gamma_{gamma}-h_{h}-eps_{eps}-dom_simple-spacing_{spacing}.bp.gz",
            "rb",
        )
    )
    ts = ts_gens(data["t"])
    v = []
    for F in data["data"]:
        v.append(
            moments.TwoLocus.Util.compute_D2_conditional(
                F, nAmin=nAmin, nAmax=nAmax, nBmin=nBmin, nBmax=nBmax
            )
            / moments.TwoLocus.Util.compute_pi2_conditional(
                F, nAmin=nAmin, nAmax=nAmax, nBmin=nBmin, nBmax=nBmax
            )
        )
    return ts, v


## Toy demography:
# t in generations, Ne ~ 1e4
nu_expand = 3.0
nu_bottle = 0.2
Ne = 1e4
t_bottle = [0, 2000, 3000, 4000]
t_expand = [0, 1000, 4000]
N_bottle = [Ne, nu_bottle * Ne, Ne]
N_expand = [Ne, nu_expand * Ne]
Ts_bottle = [0.2, 0.1, 0.1]
nus_bottle = [1, nu_bottle, 1]
Ts_expand = [0.1, 0.3]
nus_expand = [1, nu_expand]


## construct datasets
n = 50  # the sample sizes
large_thresh = -5
gammas = [-1, -10]
rhos = [0, 5]
sp = 10

## neutral data
t_neu_expand = []
sd1_neu_expand = []
sd2_neu_expand = []
data = pickle.load(
    gzip.open(
        "data/time_series/expand-n_50-rho_0.0-gamma_0-h_0.5-eps_0.0-dom_simple-spacing_10.bp.gz",
        "rb",
    )
)
for t, F in zip(data["t"], data["data"]):
    t_neu_expand.append(t)
    D = F.D()
    D2 = F.D2()
    pi2 = F.pi2()
    sd1_neu_expand.append(D / pi2)
    sd2_neu_expand.append(D2 / pi2)

tt_expand = ts_gens(data["t"])

t_neu_bottle = []
sd1_neu_bottle = []
sd2_neu_bottle = []
data = pickle.load(
    gzip.open(
        "data/time_series/bottle-n_50-rho_0.0-gamma_0-h_0.5-eps_0.0-dom_simple-spacing_10.bp.gz",
        "rb",
    )
)
for t, F in zip(data["t"], data["data"]):
    t_neu_bottle.append(t)
    D = F.D()
    D2 = F.D2()
    pi2 = F.pi2()
    sd1_neu_bottle.append(D / pi2)
    sd2_neu_bottle.append(D2 / pi2)

tt_bottle = ts_gens(data["t"])


## plots

lw = 1
colors = Category10[10]

fig = plt.figure(4, figsize=(6.5, 3.5))
fig.clf()
axes = [plt.subplot(2, 3, i) for i in range(1, 7)]

## plot the demographic models
plot_relate_curve(
    axes[0],
    t_expand,
    N_expand,
    line_style="-",
    lw=lw,
    color=c1[0],
    label="Expansion",
    gen=1,
)
plot_relate_curve(
    axes[0],
    t_bottle,
    N_bottle,
    line_style="--",
    lw=lw,
    color=c2[0],
    label="Bottleneck",
    gen=1,
)

axes[0].legend()
axes[0].set_xlim(0, 4000)
axes[0].set_xticks([0, 1000, 2000, 3000, 4000])
axes[0].set_xticklabels([4000, 3000, 2000, 1000, 0])
axes[0].set_ylim(0, 4e4)
axes[0].set_ylabel(r"$N_e$")
axes[0].set_title("Size histories")

## additivity
print("plotting additivity panel")
# gamma = -1
t, v = get_sd2(nus_expand, Ts_expand, 50, gamma=-1.0, name="expand")
axes[1].plot(t, v, "-", lw=lw, color=c1[1])
t, v = get_sd2(nus_bottle, Ts_bottle, 50, gamma=-1.0, name="bottle")
axes[1].plot(t, v, "--", lw=lw, color=c2[1])
axes[1].text(800, v[0] - 0.06, rf"$\gamma=-1$", va="center", ha="center", fontsize=7)

# gamma = -10
t, v = get_sd2(nus_expand, Ts_expand, 50, n_large=70, gamma=-10.0, name="expand")
axes[1].plot(t, v, "-", lw=lw, color=c1[1])
t, v = get_sd2(nus_bottle, Ts_bottle, 50, n_large=70, gamma=-10.0, name="bottle")
axes[1].plot(t, v, "--", lw=lw, color=c2[1])
axes[1].text(1000, v[0] + 0.04, rf"$\gamma=-10$", va="center", ha="center", fontsize=7)

# neutral baseline
axes[1].plot(tt_expand, sd2_neu_expand, "k--", lw=0.5)
axes[1].plot(tt_bottle, sd2_neu_bottle, "k--", lw=0.5)

axes[1].set_title("Interference ($\epsilon=0$)")

## synergistic epistasis
print("plotting positive epistasis panel")
# gamma = -1
t, v = get_sd2(nus_expand, Ts_expand, 50, gamma=-1.0, eps=0.5, name="expand")
axes[2].plot(t, v, "-", lw=lw, color=c1[2])
t, v = get_sd2(nus_bottle, Ts_bottle, 50, gamma=-1.0, eps=0.5, name="bottle")
axes[2].plot(t, v, "--", lw=lw, color=c2[2])
axes[2].text(800, v[0] + 0.04, rf"$\gamma=-1$", va="center", ha="center", fontsize=7)

# gamma = -10
t, v = get_sd2(
    nus_expand, Ts_expand, 50, n_large=70, gamma=-10.0, eps=0.5, name="expand"
)
axes[2].plot(t, v, "-", lw=lw, color=c1[2])
t, v = get_sd2(
    nus_bottle, Ts_bottle, 50, n_large=70, gamma=-10.0, eps=0.5, name="bottle"
)
axes[2].plot(t, v, "--", lw=lw, color=c2[2])
axes[2].text(1000, v[0] + 0.04, rf"$\gamma=-10$", va="center", ha="center", fontsize=7)

# neutral baseline
axes[2].plot(tt_expand, sd2_neu_expand, "k--", lw=0.5)
axes[2].plot(tt_bottle, sd2_neu_bottle, "k--", lw=0.5)

axes[2].set_title(r"Synergistic epistasis ($\epsilon=1/2$)")

# antagonistic epistasis
print("plotting negative epistasis panel")
# gamma = -1
t, v = get_sd2(nus_expand, Ts_expand, 50, gamma=-1.0, eps=-0.5, name="expand")
axes[3].plot(t, v, "-", lw=lw, color=c1[3])
t, v = get_sd2(nus_bottle, Ts_bottle, 50, gamma=-1.0, eps=-0.5, name="bottle")
axes[3].plot(t, v, "--", lw=lw, color=c2[3])
axes[3].text(800, v[0] - 0.06, rf"$\gamma=-1$", va="center", ha="center", fontsize=7)

# gamma = -10
t, v = get_sd2(
    nus_expand, Ts_expand, 50, n_large=70, gamma=-10.0, eps=-0.5, name="expand"
)
axes[3].plot(t, v, "-", lw=lw, color=c1[3])
t, v = get_sd2(
    nus_bottle, Ts_bottle, 50, n_large=70, gamma=-10.0, eps=-0.5, name="bottle"
)
axes[3].plot(t, v, "--", lw=lw, color=c2[3])
axes[3].text(1000, v[0] + 0.04, rf"$\gamma=-10$", va="center", ha="center", fontsize=7)

# neutral baseline
axes[3].plot(tt_expand, sd2_neu_expand, "k--", lw=0.5)
axes[3].plot(tt_bottle, sd2_neu_bottle, "k--", lw=0.5)

axes[3].set_title(r"Antagonistic epistasis ($\epsilon=-1/2$)")

# n <= 4
print("plotting rare data")
thr = 4
# each with gamma = -1
# eps = 0
t, v = get_sd2_cond(gamma=-1.0, eps=0.0, name="expand", nAmax=thr, nBmax=thr)
axes[4].plot(t, v, "-", lw=lw, color=c1[1])
t, v = get_sd2_cond(gamma=-1.0, eps=0.0, name="bottle", nAmax=thr, nBmax=thr)
axes[4].plot(t, v, "--", lw=lw, color=c2[1])
axes[4].text(2500, v[0] + 0, rf"$\epsilon=0$", va="center", ha="center", fontsize=7)
# eps = 0.5
t, v = get_sd2_cond(gamma=-1.0, eps=0.5, name="expand", nAmax=thr, nBmax=thr)
axes[4].plot(t, v, "-", lw=lw, color=c1[2])
t, v = get_sd2_cond(gamma=-1.0, eps=0.5, name="bottle", nAmax=thr, nBmax=thr)
axes[4].plot(t, v, "--", lw=lw, color=c2[2])
axes[4].text(600, v[0] - 0.03, rf"$\epsilon=1/2$", va="center", ha="center", fontsize=7)
# eps = -0.5
t, v = get_sd2_cond(gamma=-1.0, eps=-0.5, name="expand", nAmax=thr, nBmax=thr)
axes[4].plot(t, v, "-", lw=lw, color=c1[3])
t, v = get_sd2_cond(gamma=-1.0, eps=-0.5, name="bottle", nAmax=thr, nBmax=thr)
axes[4].plot(t, v, "--", lw=lw, color=c2[3])
axes[4].text(
    1000, v[0] + 0.025, rf"$\epsilon=-1/2$", va="center", ha="center", fontsize=7
)
# neutral baseline
t, v = get_sd2_cond(nAmax=thr, nBmax=thr, name="expand")
axes[4].plot(t, v, "k--", lw=0.5)
t, v = get_sd2_cond(nAmax=thr, nBmax=thr, name="bottle")
axes[4].plot(t, v, "k--", lw=0.5)

axes[4].set_title(rf"$n_A, n_B \leq {thr}$ and $\gamma=-1$")

# n >= 5
print("plotting common data")
# each with gamma = -1
# eps = 0
t, v = get_sd2_cond(gamma=-1.0, eps=0.0, name="expand", nAmin=thr + 1, nBmin=thr + 1)
axes[5].plot(t, v, "-", lw=lw, color=c1[1])
t, v = get_sd2_cond(gamma=-1.0, eps=0.0, name="bottle", nAmin=thr + 1, nBmin=thr + 1)
axes[5].plot(t, v, "--", lw=lw, color=c2[1])
axes[5].text(1000, v[0] + 0.015, rf"$\epsilon=0$", va="center", ha="center", fontsize=7)
# eps = 0.5
t, v = get_sd2_cond(gamma=-1.0, eps=0.5, name="expand", nAmin=thr + 1, nBmin=thr + 1)
axes[5].plot(t, v, "-", lw=lw, color=c1[2])
t, v = get_sd2_cond(gamma=-1.0, eps=0.5, name="bottle", nAmin=thr + 1, nBmin=thr + 1)
axes[5].plot(t, v, "--", lw=lw, color=c2[2])
axes[5].text(
    1000, v[0] - 0.025, rf"$\epsilon=1/2$", va="center", ha="center", fontsize=7
)
# eps = -0.5
t, v = get_sd2_cond(gamma=-1.0, eps=-0.5, name="expand", nAmin=thr + 1, nBmin=thr + 1)
axes[5].plot(t, v, "-", lw=lw, color=c1[3])
t, v = get_sd2_cond(gamma=-1.0, eps=-0.5, name="bottle", nAmin=thr + 1, nBmin=thr + 1)
axes[5].plot(t, v, "--", lw=lw, color=c2[3])
axes[5].text(
    1000, v[0] + 0.015, rf"$\epsilon=-1/2$", va="center", ha="center", fontsize=7
)
# neutral baseline
t, v = get_sd2_cond(nAmin=thr + 1, nBmin=thr + 1, name="expand")
axes[5].plot(t, v, "k--", lw=0.5)
t, v = get_sd2_cond(nAmin=thr + 1, nBmin=thr + 1, name="bottle")
axes[5].plot(t, v, "k--", lw=0.5)

axes[5].set_title(rf"$n_A, n_B \geq {thr + 1}$ and $\gamma=-1$")

for i in range(1, 6):
    if i > 2:
        axes[i].set_xlabel("Generations ago")
    axes[i].set_xlim(axes[0].get_xlim())
    axes[i].set_ylabel(r"$\sigma_d^2$")
    axes[i].set_xticks(axes[0].get_xticks())
    axes[i].set_xticklabels(axes[0].get_xticklabels())

fig.tight_layout()
fig.text(0.03, 0.95, "A", fontsize=8, va="center", ha="center")
fig.text(0.37, 0.95, "B", fontsize=8, va="center", ha="center")
fig.text(0.70, 0.95, "C", fontsize=8, va="center", ha="center")
fig.text(0.03, 0.49, "D", fontsize=8, va="center", ha="center")
fig.text(0.37, 0.49, "E", fontsize=8, va="center", ha="center")
fig.text(0.70, 0.49, "F", fontsize=8, va="center", ha="center")
plt.savefig("demog_bottle_expand.sd2.pdf")
