import sys, os
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

from bokeh.palettes import Colorblind


def piecewise_constant(
    nus,
    Ts,
    n,
    rho=0.0,
    gamma=0,
    h=0.5,
    dom="simple",
    eps=0.0,
    spacing=5,
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


def ts_gens(t, Ne):
    return 2 * Ne * np.array(t)


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


def get_sd1(
    nus,
    Ts,
    n,
    Ne,
    n_large=None,
    rho=0.0,
    gamma=0,
    dom="simple",
    eps=0.0,
    h=0.5,
    spacing=5,
    name=None,
):
    # see if we cached the t, sd1 vector data, so we don't need to recompute
    cached_fname = f"data/time_series/cache/{name}-n_{n}-rho_{rho}-gamma_{gamma}-h_{h}-eps_{eps}-dom_{dom}-spacing_{spacing}.sd1.bp"
    if os.path.exists(cached_fname):
        print("loading stats from cache")
        ts, v = pickle.load(open(cached_fname, "rb"))
        return ts, v
    if name is None:
        raise ValueError
    if n_large is not None:
        try:
            fname = f"data/time_series/{name}-n_{n}-rho_{rho}-gamma_{gamma}-h_{h}-eps_{eps}-dom_{dom}-spacing_{spacing}-large_{n_large}.bp.gz"
            data = pickle.load(gzip.open(fname))
        except IOError:
            print("could not find", fname)
            try:
                fname_large = f"data/time_series/{name}-n_{n_large}-rho_{rho}-gamma_{gamma}-h_{h}-eps_{eps}-dom_{dom}-spacing_{spacing}.bp.gz"
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
            nus,
            Ts,
            n,
            gamma=gamma,
            eps=eps,
            h=h,
            dom=dom,
            rho=rho,
            spacing=spacing,
            name=name,
        )
    v = []
    for F in data["data"]:
        v.append(F.D() / F.pi2())
    ts = ts_gens(data["t"], Ne)
    # cache it
    pickle.dump([ts, v], open(cached_fname, "wb+"))
    return ts, v


def get_sd1_cond(
    Ne,
    rho=0.0,
    gamma=0,
    eps=0.0,
    h=0.5,
    dom="simple",
    spacing=5,
    name=None,
    nAmin=None,
    nBmin=None,
    nAmax=None,
    nBmax=None,
):
    n = 50
    cached_fname = f"data/time_series/cache/{name}-n_{n}-rho_{rho}-gamma_{gamma}-h_{h}-eps_{eps}-dom_{dom}-spacing_{spacing}.sd1.nAmin_{nAmin}.nAmax_{nAmax}.nBmin_{nBmin}.nBmax_{nBmax}.bp"
    if os.path.exists(cached_fname):
        ts, v = pickle.load(open(cached_fname, "rb"))
        return ts, v
    data = pickle.load(
        gzip.open(
            f"data/time_series/{name}-n_{n}-rho_{rho}-gamma_{gamma}-h_{h}-eps_{eps}-dom_{dom}-spacing_{spacing}.bp.gz",
            "rb",
        )
    )
    ts = ts_gens(data["t"], Ne)
    v = []
    for F in data["data"]:
        v.append(
            moments.TwoLocus.Util.compute_D_conditional(
                F, nAmin=nAmin, nAmax=nAmax, nBmin=nBmin, nBmax=nBmax
            )
            / moments.TwoLocus.Util.compute_pi2_conditional(
                F, nAmin=nAmin, nAmax=nAmax, nBmin=nBmin, nBmax=nBmax
            )
        )
    pickle.dump([ts, v], open(cached_fname, "wb+"))
    return ts, v


## Relate demography:
pop0 = "YRI"
pop1 = "CEU"
tt, N0, N1 = load_relate_curves(
    "../../african-demography/coalescent_rate_curves/agr_refit_1_22.coal", pop0, pop1
)

# adjust to set recent size to reasonable value, set Ne for t > 1e6 years to
tt = np.concatenate(([0], tt[4:23]))
Ne = np.mean([N0[-1], N1[-1]])
N0 = N0[3:23]
N1 = N1[3:23]
Ne = np.mean([N0[-1], N1[-1]])
N0 = np.concatenate((N0, [Ne]))
N1 = np.concatenate((N1, [Ne]))

# convert to forward genetic units
Ts = ((tt[1:] - tt[:-1]) / 2 / Ne)[::-1]
nu_YRI = N0[::-1][2:] / Ne
nu_CEU = N1[::-1][2:] / Ne

## construct datasets
n = 50  # the sample sizes

## neutral data
tt_YRI, sd1_neu_YRI = get_sd1(nu_YRI, Ts, n, Ne, name="YRI")
tt_CEU, sd1_neu_CEU = get_sd1(nu_CEU, Ts, n, Ne, name="CEU")


## plots

lw = 1
colors = Colorblind[8]

fig = plt.figure(4, figsize=(6.5, 3.5))
fig.clf()
axes = [plt.subplot(2, 3, i) for i in range(1, 7)]

## plot the demographic models
gen_time = 29

plot_relate_curve(
    axes[0],
    tt * gen_time,
    N0,
    line_style="-",
    lw=lw,
    color=colors[0],
    label="YRI",
    gen=1,
)
plot_relate_curve(
    axes[0],
    tt * gen_time,
    N1,
    line_style="--",
    lw=lw,
    color=colors[1],
    label="CEU",
    gen=1,
)

axes[0].legend()
axes[0].set_xlim(1e6, 1000)
axes[0].set_xscale("log")
axes[0].set_yscale("log")
axes[0].set_xticks([1e6, 1e5, 1e4, 1e3])
axes[0].set_ylabel(r"Inverse Inst. Coal. Rate ($N_e$)")
axes[0].set_title("Size histories")

## additivity
print("plotting additivity panel")
t, v = get_sd1(nu_YRI, Ts, 50, Ne, gamma=-2.0, name="YRI")
axes[1].plot(1e6 - gen_time * t, v, "-", lw=lw, color=colors[0])
t, v = get_sd1(nu_CEU, Ts, 50, Ne, gamma=-2.0, name="CEU")
axes[1].plot(1e6 - gen_time * t, v, "--", lw=lw, color=colors[1])
# axes[1].text(2e5, v[0] + 0.02, rf"$\gamma=-1$", va="center", ha="center", fontsize=7)

# neutral baseline
axes[1].plot(1e6 - gen_time * tt_YRI, sd1_neu_YRI, "k--", lw=0.5)

axes[1].set_title("Interference")

## simple dominance
print("plotting simple dominance panel")
t, v = get_sd1(nu_YRI, Ts, 50, Ne, gamma=-2.0, h=0.1, name="YRI")
axes[2].plot(1e6 - gen_time * t, v, "-", lw=lw, color=colors[0])
t, v = get_sd1(nu_CEU, Ts, 50, Ne, gamma=-2.0, h=0.1, name="CEU")
axes[2].plot(1e6 - gen_time * t, v, "--", lw=lw, color=colors[1])
# axes[2].text(2e5, v[0] + 0.06, rf"$\gamma=-1$", va="center", ha="center", fontsize=7)

# neutral baseline
axes[2].plot(1e6 - gen_time * tt_YRI, sd1_neu_YRI, "k--", lw=0.5)

axes[2].set_title(r"Site-wise dominance ($h=0.1$)")

# gene-based dominance
print("plotting gene-based dominance panel")
t, v = get_sd1(nu_YRI, Ts, 50, Ne, gamma=-2.0, h=0.1, dom="gene", name="YRI")
axes[3].plot(1e6 - gen_time * t, v, "-", lw=lw, color=colors[0])
t, v = get_sd1(nu_CEU, Ts, 50, Ne, gamma=-2.0, h=0.1, dom="gene", name="CEU")
axes[3].plot(1e6 - gen_time * t, v, "--", lw=lw, color=colors[1])
# axes[3].text(2e5, v[0] + 0.08, rf"$\gamma=-1$", va="center", ha="center", fontsize=7)

# neutral baseline
axes[3].plot(1e6 - gen_time * tt_YRI, sd1_neu_YRI, "k--", lw=0.5)

axes[3].set_title(r"Gene-based dominance ($h=0.1$)")

# n <= 4
print("plotting rare data")
thr = 4
## each with gamma = -2
t, v = get_sd1_cond(Ne, gamma=-2.0, name="YRI", nAmax=thr, nBmax=thr)
axes[4].plot(1e6 - gen_time * t, v, "-", lw=lw, color=colors[0])
t, v = get_sd1_cond(Ne, gamma=-2.0, name="CEU", nAmax=thr, nBmax=thr)
axes[4].plot(1e6 - gen_time * t, v, "--", lw=lw, color=colors[1])
# axes[4].text(2500, v[0], rf"$\epsilon=0$", va="center", ha="center", fontsize=7)
# site-wise
t, v = get_sd1_cond(Ne, gamma=-2.0, h=0.1, name="YRI", nAmax=thr, nBmax=thr)
axes[4].plot(1e6 - gen_time * t, v, "-", lw=lw, color=colors[0])
t, v = get_sd1_cond(Ne, gamma=-2.0, h=0.1, name="CEU", nAmax=thr, nBmax=thr)
axes[4].plot(1e6 - gen_time * t, v, "--", lw=lw, color=colors[1])
# axes[4].text(2e4, v[0] - 0.5, rf"site", va="center", ha="center", fontsize=7)
# gene-based
t, v = get_sd1_cond(Ne, gamma=-2.0, h=0.1, dom="gene", name="YRI", nAmax=thr, nBmax=thr)
axes[4].plot(1e6 - gen_time * t, v, "-", lw=lw, color=colors[0])
t, v = get_sd1_cond(Ne, gamma=-2.0, h=0.1, dom="gene", name="CEU", nAmax=thr, nBmax=thr)
axes[4].plot(1e6 - gen_time * t, v, "--", lw=lw, color=colors[1])
# axes[4].text(3e5, v[0] + 1, rf"gene", va="center", ha="center", fontsize=7)
# neutral baseline
t, v = get_sd1_cond(Ne, nAmax=thr, nBmax=thr, name="YRI")
axes[4].plot(1e6 - gen_time * t, v, "k--", lw=0.5)
t, v = get_sd1_cond(Ne, nAmax=thr, nBmax=thr, name="CEU")
axes[4].plot(1e6 - gen_time * t, v, "k--", lw=0.5)

axes[4].set_title(rf"$n_A, n_B \leq {thr}$ and $\gamma=-2$")

# n >= 5
print("plotting common data")
# each with gamma = -2
t, v = get_sd1_cond(Ne, gamma=-2.0, name="YRI", nAmin=thr + 1, nBmin=thr + 1)
axes[5].plot(1e6 - gen_time * t, v, "-", lw=lw, color=colors[0])
t, v = get_sd1_cond(Ne, gamma=-2.0, name="CEU", nAmin=thr + 1, nBmin=thr + 1)
axes[5].plot(1e6 - gen_time * t, v, "--", lw=lw, color=colors[1])
# axes[5].text(2e5, v[0] - 0.2, rf"$\epsilon=0$", va="center", ha="center", fontsize=7)
# site-wise
t, v = get_sd1_cond(Ne, gamma=-2.0, h=0.1, name="YRI", nAmin=thr + 1, nBmin=thr + 1)
axes[5].plot(1e6 - gen_time * t, v, "-", lw=lw, color=colors[0])
t, v = get_sd1_cond(Ne, gamma=-2.0, h=0.1, name="CEU", nAmin=thr + 1, nBmin=thr + 1)
axes[5].plot(1e6 - gen_time * t, v, "--", lw=lw, color=colors[1])
# axes[5].text(2e5, v[0] + 0.1, rf"site", va="center", ha="center", fontsize=7)
# gene-based
t, v = get_sd1_cond(
    Ne, gamma=-2.0, h=0.1, dom="gene", name="YRI", nAmin=thr + 1, nBmin=thr + 1
)
axes[5].plot(1e6 - gen_time * t, v, "-", lw=lw, color=colors[0])
t, v = get_sd1_cond(
    Ne, gamma=-2.0, h=0.1, dom="gene", name="CEU", nAmin=thr + 1, nBmin=thr + 1
)
axes[5].plot(1e6 - gen_time * t, v, "--", lw=lw, color=colors[1])
# axes[5].text(
#    2e5, v[0] - 0.2, rf"gene", va="center", ha="center", fontsize=7
# )
# neutral baseline
t, v = get_sd1_cond(Ne, nAmin=thr + 1, nBmin=thr + 1, name="YRI")
axes[5].plot(1e6 - gen_time * t, v, "k--", lw=0.5)
t, v = get_sd1_cond(Ne, nAmin=thr + 1, nBmin=thr + 1, name="CEU")
axes[5].plot(1e6 - gen_time * t, v, "k--", lw=0.5)

axes[5].set_title(rf"$n_A, n_B \geq {thr + 1}$ and $\gamma=-1$")


for i in range(1, 6):
    if i > 2:
        axes[i].set_xlabel("Years ago")
    axes[i].set_xlim(axes[0].get_xlim())
    axes[i].set_xscale("log")
    axes[i].set_ylabel(r"$\sigma_d^1$")
    axes[i].set_xticks([1e6, 1e5, 1e4, 1e3])

fig.tight_layout()
fig.text(0.02, 0.95, "A", fontsize=8, va="center", ha="center")
fig.text(0.37, 0.95, "B", fontsize=8, va="center", ha="center")
fig.text(0.70, 0.95, "C", fontsize=8, va="center", ha="center")
fig.text(0.02, 0.49, "D", fontsize=8, va="center", ha="center")
fig.text(0.37, 0.49, "E", fontsize=8, va="center", ha="center")
fig.text(0.70, 0.49, "F", fontsize=8, va="center", ha="center")
plt.savefig("demog_YRI_CEU.dominance.pdf")
