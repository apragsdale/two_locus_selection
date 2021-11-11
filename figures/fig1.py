# Unconditioned vs allele count-conditioned LD.  Four panels, first two showing
# Hudson slices in the with rho=1 and rho=10 for neutrality and with gamma=-5
# at both loci. Second two panels show \sigma_d^2 and \sigma_d^1 for nuetrality
# and using all SNPs with gamma=-5, as well as conditioning on doubletons,
# tripletons, and few other frequencies.

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

colors = Colorblind[8]

n0 = 40  # larger sample size for improved accuracy
n = 30  # the projection size

nA = 10
nB = 10

data_fname = f"fig1_data_ns_{n0}_{n}.bp"
try:
    data1 = pickle.load(open(data_fname, "rb"))
except IOError:
    data1 = {}

if "A" not in data1.keys():
    data1["A"] = {}
    # stores the nAB distribution with n=30, nA=10, nB=10 for nuetrality
    # and gamma=-5, with rho=1 and rho=10
    for rho in [1, 10]:
        data1["A"][rho] = {}
        for gamma in [0, -5]:
            sel_params = [2 * gamma, gamma, gamma]
            F0 = moments.TwoLocus.Demographics.equilibrium(
                n0, rho=rho, sel_params=sel_params
            )
            F = F0.project(n)
            counts, pAB = moments.TwoLocus.Util.pAB(F, nA, nB)
            data1["A"][rho][gamma] = pAB
            print("computed pAB for rho, gamma =", rho, ",", gamma)
    pickle.dump(data1, open(data_fname, "wb+"))

gammas = [0, -0.1, -1, -5, -10]

if "C" not in data1.keys():
    rhos = np.logspace(-2, np.log10(30), 20)
    data1["C"] = {}
    data1["C"]["rhos"] = rhos
    data1["C"]["Ohta_Kimura"] = (5 + rhos / 2) / (11 + 13 * rhos / 2 + rhos ** 2 / 2)
    # LD curves for gamma = 0, -5, acorss rhos, for "all" SNPs, and nA=nB=2, 3, 4, 10
    for gamma in gammas:
        sel_params = [2 * gamma, gamma, gamma]
        data1["C"][gamma] = {}
        for cond in ["all", 2, 3, 4, 10]:
            data1["C"][gamma][cond] = {"sd1": [], "sd2": []}
        for rho in rhos:
            F0 = moments.TwoLocus.Demographics.equilibrium(
                n0, rho=rho, sel_params=sel_params
            )
            F = F0.project(n)
            for cond in data["C"][gamma].keys():
                if cond == "all":
                    data1["C"][gamma]["all"]["sd1"].append(F.D() / F.pi2())
                    data1["C"][gamma]["all"]["sd2"].append(F.D2() / F.pi2())
                else:
                    data1["C"][gamma][cond]["sd1"].append(
                        F.D(nA=cond, nB=cond) / F.pi2(nA=cond, nB=cond)
                    )
                    data1["C"][gamma][cond]["sd2"].append(
                        F.D2(nA=cond, nB=cond) / F.pi2(nA=cond, nB=cond)
                    )
            print("finished gamma, rho = ", gamma, ",", rho)
    pickle.dump(data1, open(data_fname, "wb+"))

markers = [".", "x", "+", "1", "2"]
rhos = data1["C"]["rhos"]

ms = 4

fig = plt.figure(1, figsize=(6.5, 5))

## previous c and d from figure 1

ax1 = plt.subplot(3, 2, 1)

ax1.plot(rhos, data1["C"]["Ohta_Kimura"], "k--", lw=1, label="Ohta-Kimura")

for ii, gamma in enumerate(gammas):
    ax1.plot(
        rhos,
        data1["C"][gamma]["all"]["sd2"],
        markers[ii] + "-",
        ms=ms,
        lw=0.5,
        label=rf"$\gamma={gamma}$",
        color=colors[ii],
    )

ax1.set_yscale("log")
ax1.set_xscale("log")
ax1.set_ylabel(r"$\sigma_d^2$")
ax1.set_title("Normalized variance of $D$ ($\sigma_d^2$)")

ax2 = plt.subplot(3, 2, 2)

ax2.plot(rhos, 0 * data1["C"]["Ohta_Kimura"], "k--", lw=1, label=None)

for ii, gamma in enumerate(gammas):
    ax2.plot(
        rhos,
        data1["C"][gamma]["all"]["sd1"],
        markers[ii] + "-",
        ms=ms,
        lw=0.5,
        label=rf"$\gamma={gamma}$",
        color=colors[ii],
    )

ax2.set_xscale("log")
ax2.set_ylabel(r"$\sigma_d^1$")
ax2.set_title("Normalized signed LD ($\sigma_d^1$)")
ax2.legend()


### previous figure 2
epsilons = [-0.5, 0.0, 0.5, 1.0]
gammas = [-1, -10]
conditions = ["all", 2, 3, 4, 10]

n0 = 50

data_fname = f"fig1b_data_ns_{n0}_{n}.bp"
try:
    data2 = pickle.load(open(data_fname, "rb"))
except IOError:
    data2 = {}
    # stores the sigma_d^2 and sigma_d^1 curves for all snps and
    rhos = np.logspace(-2, np.log10(30), 20)
    data2["rhos"] = rhos
    data2["Ohta_Kimura"] = (5 + rhos / 2) / (11 + 13 * rhos / 2 + rhos ** 2 / 2)
    for gamma in gammas:
        for eps in epsilons:
            print("Running with  gamma, eps = ", gamma, ",", eps, end=" ")
            sel_params = moments.TwoLocus.Util.additive_epistasis(gamma, epsilon=eps)
            data2[(gamma, eps)] = {cond: {"sd1": [], "sd2": []} for cond in conditions}
            for rho in rhos:
                print(".", end="")
                F0 = moments.TwoLocus.Demographics.equilibrium(
                    n0, rho=rho, sel_params=sel_params
                )
                F = F0.project(n)
                for cond in conditions:
                    if cond == "all":
                        data2[(gamma, eps)][cond]["sd1"].append(F.D() / F.pi2())
                        data2[(gamma, eps)][cond]["sd2"].append(F.D2() / F.pi2())
                    else:
                        data2[(gamma, eps)][cond]["sd1"].append(
                            F.D(nA=cond, nB=cond) / F.pi2(nA=cond, nB=cond)
                        )
                        data2[(gamma, eps)][cond]["sd2"].append(
                            F.D2(nA=cond, nB=cond) / F.pi2(nA=cond, nB=cond)
                        )
            print("finished")
            print("sd1:", data2[(gamma, eps)]["all"]["sd1"])
            print("sd2:", data2[(gamma, eps)]["all"]["sd2"])
            print()
    pickle.dump(data2, open(data_fname, "wb+"))


rhos = data2["rhos"]
ok = data2["Ohta_Kimura"]

markers = ["x", "+", "1", "2"]
colors = Colorblind[8]

# top panels are sigma_d^2
ax3 = plt.subplot(3, 2, 3)
ax4 = plt.subplot(3, 2, 4)

ax3.plot(rhos, ok, "k--", lw=1, label=None)
ax4.plot(rhos, 0 * ok, "k--", lw=1, label=None)

for ii, eps in enumerate(epsilons):
    ax3.plot(
        rhos,
        data2[(-1, eps)]["all"]["sd2"],
        markers[ii] + "-",
        color=colors[ii],
        ms=4,
        lw=0.5,
        label=rf"$\epsilon={eps}$",
    )

for ii, eps in enumerate(epsilons):
    ax4.plot(
        rhos,
        [2 * x for x in data2[(-1, eps)]["all"]["sd1"]],
        markers[ii] + "-",
        color=colors[ii],
        ms=4,
        lw=0.5,
        label=rf"$\epsilon={eps}$",
    )


ax3.set_xscale("log")
ax4.set_xscale("log")
ax3.set_yscale("log")
ax3.set_ylabel(r"$\sigma_d^2$")
ax4.set_ylabel(r"$\sigma_d^1$")
ax3.legend()
ax3.set_title(r"$\gamma = -1$")
ax4.set_title(r"$\gamma = -1$")

ax5 = plt.subplot(3, 2, 5)
ax6 = plt.subplot(3, 2, 6)

ax5.plot(rhos, ok, "k--", lw=1, label=None)
ax6.plot(rhos, 0 * ok, "k--", lw=1, label=None)

for ii, eps in enumerate(epsilons):
    ax5.plot(
        rhos,
        data2[(-10, eps)]["all"]["sd2"],
        markers[ii] + "-",
        color=colors[ii],
        ms=4,
        lw=0.5,
        label=rf"$\epsilon={eps}$",
    )


for ii, eps in enumerate(epsilons):
    ax6.plot(
        rhos,
        [2 * x for x in data2[(-10, eps)]["all"]["sd1"]],
        markers[ii] + "-",
        color=colors[ii],
        ms=4,
        lw=0.5,
        label=rf"$\epsilon={eps}$",
    )

ax5.set_xscale("log")
ax6.set_xscale("log")
ax5.set_yscale("log")
ax5.set_ylabel(r"$\sigma_d^2$")
ax6.set_ylabel(r"$\sigma_d^1$")
ax5.legend()
ax5.set_title(r"$\gamma = -10$")
ax6.set_title(r"$\gamma = -10$")
ax5.set_xlabel(r"$\rho$")
ax6.set_xlabel(r"$\rho$")

ax3.set_ylim([.005, 0.5])
ax4.set_ylim([-1.2, 1.2])
ax5.set_ylim(ax3.get_ylim())
ax6.set_ylim(ax4.get_ylim())

fig.tight_layout()
fig.text(0.02, 0.96, "A", fontsize=8, ha="center", va="center")
fig.text(0.52, 0.96, "B", fontsize=8, ha="center", va="center")
fig.text(0.02, 0.65, "C", fontsize=8, ha="center", va="center")
fig.text(0.52, 0.65, "D", fontsize=8, ha="center", va="center")
fig.text(0.02, 0.32, "E", fontsize=8, ha="center", va="center")
fig.text(0.52, 0.32, "F", fontsize=8, ha="center", va="center")

plt.savefig(f"fig1_combo.pdf")
plt.savefig(f"fig1_combo.png", dpi=300)
