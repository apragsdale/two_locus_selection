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
    data = pickle.load(open(data_fname, "rb"))
except IOError:
    data = {}

if "A" not in data.keys():
    data["A"] = {}
    # stores the nAB distribution with n=30, nA=10, nB=10 for nuetrality
    # and gamma=-5, with rho=1 and rho=10
    for rho in [1, 10]:
        data["A"][rho] = {}
        for gamma in [0, -5]:
            sel_params = [2 * gamma, gamma, gamma]
            F0 = moments.TwoLocus.Demographics.equilibrium(
                n0, rho=rho, sel_params=sel_params
            )
            F = F0.project(n)
            counts, pAB = moments.TwoLocus.Util.pAB(F, nA, nB)
            data["A"][rho][gamma] = pAB
            print("computed pAB for rho, gamma =", rho, ",", gamma)
    pickle.dump(data, open(data_fname, "wb+"))

gammas = [0, -0.1, -1, -5, -10]

if "C" not in data.keys():
    rhos = np.logspace(-2, np.log10(30), 20)
    data["C"] = {}
    data["C"]["rhos"] = rhos
    data["C"]["Ohta_Kimura"] = (5 + rhos / 2) / (11 + 13 * rhos / 2 + rhos ** 2 / 2)
    # LD curves for gamma = 0, -5, acorss rhos, for "all" SNPs, and nA=nB=2, 3, 4, 10
    for gamma in gammas:
        sel_params = [2 * gamma, gamma, gamma]
        data["C"][gamma] = {}
        for cond in ["all", 2, 3, 4, 10]:
            data["C"][gamma][cond] = {"sd1": [], "sd2": []}
        for rho in rhos:
            F0 = moments.TwoLocus.Demographics.equilibrium(
                n0, rho=rho, sel_params=sel_params
            )
            F = F0.project(n)
            for cond in data["C"][gamma].keys():
                if cond == "all":
                    data["C"][gamma]["all"]["sd1"].append(F.D() / F.pi2())
                    data["C"][gamma]["all"]["sd2"].append(F.D2() / F.pi2())
                else:
                    data["C"][gamma][cond]["sd1"].append(
                        F.D(nA=cond, nB=cond) / F.pi2(nA=cond, nB=cond)
                    )
                    data["C"][gamma][cond]["sd2"].append(
                        F.D2(nA=cond, nB=cond) / F.pi2(nA=cond, nB=cond)
                    )
            print("finished gamma, rho = ", gamma, ",", rho)
    pickle.dump(data, open(data_fname, "wb+"))

fig = plt.figure(1, figsize=(6.5, 4))
fig.clf()

ax1 = plt.subplot(2, 2, 1)

ax1.bar(
    np.arange(0, 11) - 0.2,
    data["A"][1][0] / data["A"][1][0].sum(),
    width=0.35,
    label=r"$\gamma=0$",
    color=colors[0],
)
ax1.bar(
    np.arange(0, 11) + 0.2,
    data["A"][1][-5] / data["A"][1][-5].sum(),
    width=0.35,
    label=r"$\gamma=-5$",
    color=colors[1],
)
ax1.set_title(rf"$\rho = 1$, $n=30$, $n_A=n_B=10$")
ax1.set_ylabel("Probability")
ax1.set_xlabel(r"$n_{AB}$")
ax1.legend(frameon=False)

ax2 = plt.subplot(2, 2, 2)

ax2.bar(
    np.arange(0, 11) - 0.2,
    data["A"][10][0] / data["A"][10][0].sum(),
    width=0.35,
    label=r"$\gamma=0$",
    color=colors[0],
)
ax2.bar(
    np.arange(0, 11) + 0.2,
    data["A"][10][-5] / data["A"][10][-5].sum(),
    width=0.35,
    label=r"$\gamma=-5$",
    color=colors[1],
)
ax2.set_title(rf"$\rho = 10$, $n=30$, $n_A=n_B=10$")
ax2.set_ylabel("Probability")
ax2.set_xlabel(r"$n_{AB}$")

markers = [".", "x", "+", "1", "2"]
rhos = data["C"]["rhos"]

ms = 4

ax3 = plt.subplot(2, 2, 3)

ax3.plot(rhos, data["C"]["Ohta_Kimura"], "k--", lw=1, label="Ohta-Kimura")

for ii, gamma in enumerate(gammas):
    ax3.plot(
        rhos,
        data["C"][gamma]["all"]["sd2"],
        markers[ii] + "-",
        ms=ms,
        lw=0.5,
        label=rf"$\gamma={gamma}$",
        color=colors[ii],
    )

ax3.set_yscale("log")
ax3.set_xscale("log")
ax3.set_ylabel(r"$\sigma_d^2$")
ax3.set_xlabel(r"$\rho$")
#ax3.legend(ncol=2)

ax4 = plt.subplot(2, 2, 4)

ax4.plot(rhos, 0 * data["C"]["Ohta_Kimura"], "k--", lw=1, label=None)

for ii, gamma in enumerate(gammas):
    ax4.plot(
        rhos,
        data["C"][gamma]["all"]["sd1"],
        markers[ii] + "-",
        ms=ms,
        lw=0.5,
        label=rf"$\gamma={gamma}$",
        color=colors[ii],
    )

ax4.set_xscale("log")
ax4.set_ylabel(r"$\sigma_d^1$")
ax4.set_xlabel(r"$\rho$")
ax4.legend()

fig.tight_layout()
fig.text(0.02, 0.98, "A", fontsize=8, ha="center", va="center")
fig.text(0.52, 0.98, "B", fontsize=8, ha="center", va="center")
fig.text(0.02, 0.48, "C", fontsize=8, ha="center", va="center")
fig.text(0.52, 0.48, "D", fontsize=8, ha="center", va="center")

plt.savefig(f"fig1_n0_{n0}.pdf")
