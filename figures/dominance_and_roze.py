# Unconditioned vs allele count-conditioned LD.  Four panels, first two showing
# Hudson slices in the with rho=1 and rho=10 for neutrality and with gamma=-5
# at both loci. Second two panels show \sigma_d^2 and \sigma_d^1 for nuetrality
# and using all SNPs with gamma=-5, as well as conditioning on doubletons,
# tripletons, and few other frequencies.

import moments.TwoLocus
import matplotlib.pylab as plt
import numpy as np
import pickle, gzip

# set font sizes
import matplotlib

plt.rcParams["legend.title_fontsize"] = "xx-small"
matplotlib.rc("xtick", labelsize=6)
matplotlib.rc("ytick", labelsize=6)
matplotlib.rc("axes", labelsize=7)
matplotlib.rc("axes", titlesize=7)
matplotlib.rc("legend", fontsize=6)

from bokeh.palettes import Colorblind

n0 = 70  # larger sample size for improved accuracy
n = 30  # the projection size

hs = [0.0, 0.1, 0.5, 1.0]
gammas = [-1, -5]
conditions = ["all", 2, 3, 4, 10]

data_fname = f"data/fig2_data_ns_{n0}_{n}.bp"
try:
    data = pickle.load(open(data_fname, "rb"))
except IOError:
    data = {}
    data["ns"] = (n0, n)
    # stores the sigma_d^2 and sigma_d^1 curves for all snps and
    rhos = np.logspace(-2, np.log10(30), 20)
    data["rhos"] = rhos
    data["Ohta_Kimura"] = (5 + rhos / 2) / (11 + 13 * rhos / 2 + rhos ** 2 / 2)
    for gamma in gammas:
        for h in hs:
            print("Running gamma, h = ", gamma, ",", h, end=" ", flush=True)
            sel_params = moments.TwoLocus.Util.simple_dominance(gamma, h=h)
            data[(gamma, h)] = {cond: {"sd1": [], "sd2": []} for cond in conditions}
            for rho in rhos:
                print(".", end="", flush=True)
                F0 = moments.TwoLocus.Demographics.equilibrium(
                    n0, rho=rho, sel_params_general=sel_params
                )
                F = F0.project(n)
                for cond in conditions:
                    if cond == "all":
                        data[(gamma, h)][cond]["sd1"].append(F.D() / F.pi2())
                        data[(gamma, h)][cond]["sd2"].append(F.D2() / F.pi2())
                    else:
                        data[(gamma, h)][cond]["sd1"].append(
                            F.D(nA=cond, nB=cond) / F.pi2(nA=cond, nB=cond)
                        )
                        data[(gamma, h)][cond]["sd2"].append(
                            F.D2(nA=cond, nB=cond) / F.pi2(nA=cond, nB=cond)
                        )
            print("finished.")
            print("sd1:", data[(gamma, h)]["all"]["sd1"])
            print("sd2:", data[(gamma, h)]["all"]["sd2"])
            print(flush=True)
    pickle.dump(data, open(data_fname, "wb+"))


rhos = data["rhos"]
ok = data["Ohta_Kimura"]

colors = Colorblind[8]

fig = plt.figure(19382, figsize=(6.5, 4.5))
fig.clf()

grid = (5, 6)
# sigma_d^1 for gamma = -1 and -5
ax1 = plt.subplot2grid(grid, (0, 0), rowspan=2, colspan=2)
ax2 = plt.subplot2grid(grid, (0, 2), rowspan=2, colspan=2)
# roze comparisons
ax3 = plt.subplot2grid(grid, (2, 0), rowspan=3, colspan=3)
ax4 = plt.subplot2grid(grid, (2, 3), rowspan=3, colspan=3)
# varying gamma (upper right)
ax5 = plt.subplot2grid(grid, (0, 4), rowspan=2, colspan=2)

ax1.plot(rhos, 0 * ok, "k--", lw=1, label=None)
ax2.plot(rhos, 0 * ok, "k--", lw=1, label=None)

for ii, h in enumerate(hs):
    ax1.plot(
        rhos,
        [2 * x for x in data[(-1, h)]["all"]["sd1"]],
        color=colors[ii],
        lw=1,
        label=rf"$h={h}$",
    )

for ii, h in enumerate(hs):
    ax2.plot(
        rhos,
        [2 * x for x in data[(-5, h)]["all"]["sd1"]],
        color=colors[ii],
        lw=1,
        label=rf"$h={h}$",
    )


data_vary = pickle.load(gzip.open("data/dominance_vary_gamma.bp.gz", "rb"))
gammas = data_vary["gammas"]
hs = sorted(data_vary["sd1"].keys())

ax5.plot(-gammas, 0*gammas, "k--", lw=1, label=None)
for ii, h in enumerate(hs):
    ax5.plot(-gammas, data_vary["sd1"][h], color=colors[ii + 4], lw=1, label=rf"$h={h}$")

ax5.set_title(r"$\rho=0$")
ax5.set_xlabel(r"$\gamma=2Ns$")
ax5.set_xscale("log")
ax5.legend(fontsize=6, loc="lower left")
ax5.set_xticks([0.1, 1, 10])
ax5.set_xticklabels([-0.1, -1, -10])

# roze data
N = 1000
n = 50

hs = np.linspace(0.05, 0.5, 10)
rs = [0.001, 0.0001]
shs = [0.01, 0.001, 0.0004, 0.0001]

vals = pickle.load(open("data/roze_sd1s.bp", "rb"))
exp = pickle.load(open("data/roze_comp_data.bp", "rb"))

r = 0.0001
for sh in shs[::-1]:
    v = [vals[r][sh][h] for h in hs]
    ax3.plot(hs, v, "o--", ms=3, lw=0.5, label=f"$-{sh}$ (sim.)")

ax3.set_prop_cycle(None)
for sh in shs[::-1]:
    v = [exp[4 * N * r][sh][h] for h in sorted(exp[4 * N * r][sh].keys())]
    if sh == 0.01:
        continue
    elif sh == 0.001:
        ax3.plot(hs[1:], v[1:], "x-", ms=5, lw=1, label="moments")
    else:
        ax3.plot(hs, v, "x-", ms=5, lw=1, label="moments")

r = 0.001
for sh in shs[::-1]:
    v = [vals[r][sh][h] for h in hs]
    ax4.plot(hs, v, "o--", ms=3, lw=0.5, label=f"$-{sh}$ (sim.)")

ax4.set_prop_cycle(None)
for sh in shs[::-1]:
    v = [exp[4 * N * r][sh][h] for h in sorted(exp[4 * N * r][sh].keys())]
    if sh == 0.01:
        continue
    elif sh == 0.001:
        ax4.plot(hs[1:], v[1:], "x-", ms=5, lw=1, label="moments")
    else:
        ax4.plot(hs, v, "x-", ms=5, lw=1, label="moments")


ax1.set_xscale("log")
ax2.set_xscale("log")
ax1.set_ylabel(r"$\sigma_d^1$")
#ax2.set_ylabel(r"$\sigma_d^1$")
ax1.legend()
ax2.legend()
ax1.set_title(r"$2Ns = -1$")
ax2.set_title(r"$2Ns = -5$")
ax1.set_xlabel(r"$\rho$")
ax2.set_xlabel(r"$\rho$")
ax1.set_ylim([-0.4, 0.1])
ax2.set_ylim([-1.2, 0.1])

ax3.set_ylabel(r"$\sigma_d^1$")
#ax4.set_ylabel(r"$\sigma_d^1$")
ax3.legend(title="$sh$", ncol=2, fontsize=6, loc="lower right")
ax3.set_xlim(0, 0.55)
ax3.set_ylim(-0.51, 0.11)
ax4.legend(title="$sh$", ncol=2, fontsize=6, loc="lower right")
ax4.set_xlim(0, 0.55)
ax4.set_ylim(-0.11, 0.05)
ax3.set_title(r"$r = 0.0001, N=1000$")
ax4.set_title(r"$r = 0.001, N=1000$")
ax3.set_xlabel(r"$h$")
ax4.set_xlabel(r"$h$")

fig.tight_layout()
fig.text(0.04, 0.97, "A", fontsize=8, ha="center", va="center")
fig.text(0.37, 0.97, "B", fontsize=8, ha="center", va="center")
fig.text(0.68, 0.97, "C", fontsize=8, ha="center", va="center")
fig.text(0.04, 0.57, "D", fontsize=8, ha="center", va="center")
fig.text(0.53, 0.57, "E", fontsize=8, ha="center", va="center")

plt.savefig(f"dominance_roze.pdf")
plt.savefig(f"dominance_roze.png", dpi=300)
