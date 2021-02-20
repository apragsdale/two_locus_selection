# Figure 3: Demography and trajectories.
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

try:
    n = int(sys.argv[1])  # the sample size
    gamma = float(sys.argv[2])
    rho = float(sys.argv[3])
except:
    raise ValueError("inputs are n, gamma, rho in that order")

def piecewise_constant(
    nus,
    Ts,
    n,
    rho=0.0,
    gamma=0,
    h=None,
    eps=None,
    track_stats=True,
    spacing=3,
    conditions=[],
):
    if track_stats:
        stats = {"t": [], "sd1": {"all": [],}, "sd2": {"all": [],}}
        for cond in conditions:
            stats["sd1"][cond] = []
            stats["sd2"][cond] = []
    if gamma == 0:
        sel_params = None
        sel_params_general = None
    elif h is None and eps is None:
        sel_params = moments.TwoLocus.Util.additive_epistasis(gamma)
        sel_params_general = None
    elif h is not None:
        sel_params = None
        sel_params_general = moments.TwoLocus.Util.simple_dominance(gamma, h=h)
    elif eps is not None:
        sel_params = moments.TwoLocus.Util.additive_epistasis(gamma, epsilon=eps)
        sel_params_general = None
    F = moments.TwoLocus.Demographics.equilibrium(
        n, rho=rho, sel_params=sel_params, sel_params_general=sel_params_general
    )
    if track_stats:
        stats["t"].append(0)
        stats["sd1"]["all"].append(F.D() / F.pi2())
        stats["sd2"]["all"].append(F.D2() / F.pi2())
        for cond in conditions:
            stats["sd1"][cond].append(F.D(nA=cond, nB=cond) / F.pi2(nA=cond, nB=cond))
            stats["sd2"][cond].append(F.D2(nA=cond, nB=cond) / F.pi2(nA=cond, nB=cond))
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
                stats["t"].append(stats["t"][-1] + T_sub)
                stats["sd1"]["all"].append(F.D() / F.pi2())
                stats["sd2"]["all"].append(F.D2() / F.pi2())
                for cond in conditions:
                    stats["sd1"][cond].append(
                        F.D(nA=cond, nB=cond) / F.pi2(nA=cond, nB=cond)
                    )
                    stats["sd2"][cond].append(
                        F.D2(nA=cond, nB=cond) / F.pi2(nA=cond, nB=cond)
                    )
    print()
    return F, stats


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


data = {}
data[pop0] = {}
data[pop1] = {}
# keys are (rho, gamma, eps, h, gene-based) where gene-based is True or False

print(f"computing {pop0}, {pop1} stats for rho={rho}, gamma=0")
F0, stats0 = piecewise_constant(nu0s, Ts, n, rho=rho, gamma=0, conditions=[2])
F1, stats1 = piecewise_constant(nu1s, Ts, n, rho=rho, gamma=0, conditions=[2])
data[pop0][(rho, 0, 0, 0.5, False)] = stats0
data[pop1][(rho, 0, 0, 0.5, False)] = stats1

eps = 0.5
print(f"computing {pop0}, {pop1} stats for rho={rho}, gamma={gamma}, eps={eps}")
F0, stats0 = piecewise_constant(
    nu0s, Ts, n, rho=rho, gamma=gamma, eps=eps, conditions=[2]
)
F1, stats1 = piecewise_constant(
    nu1s, Ts, n, =rho=rho, gamma=gamma, eps=eps, conditions=[2]
)
data[pop0][(rho, gamma, eps, 0.5, False)] = stats0
data[pop1][(rho, gamma, eps, 0.5, False)] = stats1

print(f"computing {pop0}, {pop1} stats for rho={rho}, gamma={gamma}")
F0, stats0 = piecewise_constant(nu0s, Ts, n, rho=rho, gamma=gamma, conditions=[2])
F1, stats1 = piecewise_constant(nu1s, Ts, n, rho=rho, gamma=gamma, conditions=[2])
data[pop0][(rho, gamma, 0, 0.5, False)] = stats0
data[pop1][(rho, gamma, 0, 0.5, False)] = stats1

eps = 0.5
print(f"computing {pop0}, {pop1} stats for rho={rho}, gamma={gamma}, eps={eps}")
F0, stats0 = piecewise_constant(
    nu0s, Ts, n, rho=rho, gamma=gamma, eps=eps, conditions=[2]
)
F1, stats1 = piecewise_constant(
    nu1s, Ts, n, rho=rho, gamma=gamma, eps=eps, conditions=[2]
)
data[pop0][(rho, gamma, eps, 0.5, False)] = stats0
data[pop1][(rho, gamma, eps, 0.5, False)] = stats1


h = 0.1
print(f"computing {pop0}, {pop1} stats for rho={rho}, gamma={gamma}, h={h}")
F0, stats0 = piecewise_constant(nu0s, Ts, n, rho=rho, gamma=gamma, h=h, conditions=[2])
F1, stats1 = piecewise_constant(nu1s, Ts, n, rho=rho, gamma=gamma, h=h, conditions=[2])
data[pop0][(rho, gamma, 0, h, False)] = stats0
data[pop1][(rho, gamma, 0, h, False)] = stats1

Ts_years = 1e6 - 2 * Ne * 29 * np.array(stats0["t"])

fig = plt.figure(4, figsize=(3.25, 6))
fig.clf()

markers = ["x", "+", "1"]
colors = Colorblind[8]

ax1 = plt.subplot(3, 1, 1)

plot_relate_curve(ax1, t, N0, line_style="-", lw=1, color=colors[0], label=pop0)
plot_relate_curve(ax1, t, N1, line_style="--", lw=1, color=colors[1], label=pop1)

ax1.legend(frameon=False)
ax1.set_xscale("log")
ax1.set_yscale("log")
ax1.set_xlim(1e6, 1e3)
ax1.set_ylim(1e2, 1e7)
ax1.set_ylabel(r"$N_e$")

ax2 = plt.subplot(3, 1, 2)
ax3 = plt.subplot(3, 1, 3)

ax2.plot(
    Ts_years,
    data[pop0][(rho, 0, 0, 0.5, False)]["sd1"]["all"],
    label=rf"{pop0}, $\gamma=0$",
    color="black",
    lw=0.5,
)
ax2.plot(
    Ts_years,
    data[pop1][(rho, 0, 0, 0.5, False)]["sd1"]["all"],
    label=rf"{pop1}, $\gamma=0$",
    linestyle="--",
    color="black",
    lw=0.5,
)
ax2.plot(
    Ts_years,
    data[pop0][(rho, gamma, 0, 0.5, False)]["sd1"]["all"],
    label=rf"{pop0}, $\gamma={gamma}$",
    color=colors[2],
    lw=1,
)
ax2.plot(
    Ts_years,
    data[pop1][(rho, gamma, 0, 0.5, False)]["sd1"]["all"],
    label=rf"{pop1}, $\gamma={gamma}$",
    linestyle="--",
    color=colors[2],
    lw=1,
)
ax2.plot(
    Ts_years,
    data[pop0][(rho, gamma, 0.5, 0.5, False)]["sd1"]["all"],
    label=rf"{pop0}, $\epsilon={eps}$",
    color=colors[3],
    lw=1,
)
ax2.plot(
    Ts_years,
    data[pop1][(rho, gamma, 0.5, 0.5, False)]["sd1"]["all"],
    label=rf"{pop1}, $\epsilon={eps}$",
    linestyle="--",
    color=colors[3],
    lw=1,
)
ax2.plot(
    Ts_years,
    data[pop0][(rho, gamma, 0, 0.1, False)]["sd1"]["all"],
    label=rf"{pop0}, $h={h}$",
    color=colors[4],
    lw=1,
)
ax2.plot(
    Ts_years,
    data[pop1][(rho, gamma, 0, 0.1, False)]["sd1"]["all"],
    label=rf"{pop1}, $h={h}$",
    linestyle="--",
    color=colors[4],
    lw=1,
)


ax3.plot(
    Ts_years,
    data[pop0][(rho, 0, 0, 0.5, False)]["sd1"][2],
    label=rf"{pop0}, $\gamma=0$",
    color="black",
    lw=0.5,
)
ax3.plot(
    Ts_years,
    data[pop1][(rho, 0, 0, 0.5, False)]["sd1"][2],
    label=rf"{pop1}, $\gamma=0$",
    linestyle="--",
    color="black",
    lw=0.5,
)
ax3.plot(
    Ts_years,
    data[pop0][(rho, gamma, 0, 0.5, False)]["sd1"][2],
    label=rf"{pop0}, $\gamma={gamma}$",
    color=colors[2],
    lw=1,
)
ax3.plot(
    Ts_years,
    data[pop1][(rho, gamma, 0, 0.5, False)]["sd1"][2],
    label=rf"{pop1}, $\gamma={gamma}$",
    linestyle="--",
    color=colors[2],
    lw=1,
)
ax3.plot(
    Ts_years,
    data[pop0][(rho, gamma, 0.5, 0.5, False)]["sd1"][2],
    label=rf"{pop0}, $\epsilon={eps}$",
    color=colors[3],
    lw=1,
)
ax3.plot(
    Ts_years,
    data[pop1][(rho, gamma, 0.5, 0.5, False)]["sd1"][2],
    label=rf"{pop1}, $\epsilon={eps}$",
    linestyle="--",
    color=colors[3],
    lw=1,
)
ax3.plot(
    Ts_years,
    data[pop0][(rho, gamma, 0, 0.1, False)]["sd1"][2],
    label=rf"{pop0}, $h={h}$",
    color=colors[4],
    lw=1,
)
ax3.plot(
    Ts_years,
    data[pop1][(rho, gamma, 0, 0.1, False)]["sd1"][2],
    label=rf"{pop1}, $h={h}$",
    linestyle="--",
    color=colors[4],
    lw=1,
)

ax2.set_xlim(ax1.get_xlim())
ax2.set_xscale("log")
ax3.set_xlim(ax1.get_xlim())
ax3.set_xscale("log")

ax2.legend(frameon=False, fontsize=5, ncol=2)
ax3.legend(frameon=False, fontsize=5, ncol=2)

ax2.set_ylabel(f"$\sigma_d^1$")
ax3.set_ylabel(f"$\sigma_d^1$")

ax3.set_xlabel("Years ago")
ax2.set_title("All frequencies")
ax3.set_title("$n_A = n_B = 2$")

fig.tight_layout()

# fig.text(0.05, 0.97, "A", fontsize=8, ha="center", va="center")

plt.savefig(f"fig4_n0_{n}_rho_{rho}_gamma_{gamma}.pdf")
