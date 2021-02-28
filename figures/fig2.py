# Figure 2: Additive selection with epistasis.  In this epistasis model, sAB =
# (1 + epsilon) * (sA + sB).  Smaller, single-column figure showing the effect
# of synergistic epistasis for weakly and moderately selected variants.
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

## got back and run with larger n0
n0 = int(sys.argv[1])  # larger sample size for improved accuracy
n = 30  # the projection size

epsilons = [-0.5, 0.0, 0.5, 1.0]
gammas = [-1, -10]
conditions = ["all", 2, 3, 4, 10]

data_fname = f"fig2_data_ns_{n0}_{n}.bp"
try:
    data = pickle.load(open(data_fname, "rb"))
except IOError:
    data = {}
    # stores the sigma_d^2 and sigma_d^1 curves for all snps and
    rhos = np.logspace(-2, np.log10(30), 20)
    data["rhos"] = rhos
    data["Ohta_Kimura"] = (5 + rhos / 2) / (11 + 13 * rhos / 2 + rhos ** 2 / 2)
    for gamma in gammas:
        for eps in epsilons:
            print("Running with  gamma, eps = ", gamma, ",", eps, end=" ")
            sel_params = moments.TwoLocus.Util.additive_epistasis(gamma, epsilon=eps)
            data[(gamma, eps)] = {cond: {"sd1": [], "sd2": []} for cond in conditions}
            for rho in rhos:
                print(".", end="")
                F0 = moments.TwoLocus.Demographics.equilibrium(
                    n0, rho=rho, sel_params=sel_params
                )
                F = F0.project(n)
                for cond in conditions:
                    if cond == "all":
                        data[(gamma, eps)][cond]["sd1"].append(F.D() / F.pi2())
                        data[(gamma, eps)][cond]["sd2"].append(F.D2() / F.pi2())
                    else:
                        data[(gamma, eps)][cond]["sd1"].append(
                            F.D(nA=cond, nB=cond) / F.pi2(nA=cond, nB=cond)
                        )
                        data[(gamma, eps)][cond]["sd2"].append(
                            F.D2(nA=cond, nB=cond) / F.pi2(nA=cond, nB=cond)
                        )
            print("finished")
            print("sd1:", data[(gamma, eps)]["all"]["sd1"])
            print("sd2:", data[(gamma, eps)]["all"]["sd2"])
            print()
    pickle.dump(data, open(data_fname, "wb+"))


fig = plt.figure(2, figsize=(6.5, 3.5))
fig.clf()

rhos = data["rhos"]
ok = data["Ohta_Kimura"]

markers = ["x", "+", "1", "2"]
colors = Colorblind[8]

# top panels are sigma_d^2
ax1 = plt.subplot(2, 2, 1)
ax2 = plt.subplot(2, 2, 2)

ax1.plot(rhos, ok, "k--", lw=1, label=None)
ax2.plot(rhos, ok, "k--", lw=1, label=None)

for ii, eps in enumerate(epsilons):
    ax1.plot(
        rhos,
        data[(-1, eps)]["all"]["sd2"],
        markers[ii] + "-",
        color=colors[ii],
        ms=4,
        lw=0.5,
        label=rf"$\epsilon={eps}$",
    )

for ii, eps in enumerate(epsilons):
    ax2.plot(
        rhos,
        data[(-10, eps)]["all"]["sd2"],
        markers[ii] + "-",
        color=colors[ii],
        ms=4,
        lw=0.5,
        label=rf"$\epsilon={eps}$",
    )

ax1.set_xscale("log")
ax1.set_yscale("log")
ax2.set_xscale("log")
ax2.set_yscale("log")
ax1.set_ylabel(r"$\sigma_d^2$")
ax1.legend()
#ax2.legend()
ax1.set_title(r"$\gamma = -1$")
ax2.set_title(r"$\gamma = -10$")

ax3 = plt.subplot(2, 2, 3)
ax4 = plt.subplot(2, 2, 4)

ax3.plot(rhos, 0 * ok, "k--", lw=1, label=None)
ax4.plot(rhos, 0 * ok, "k--", lw=1, label=None)

for ii, eps in enumerate(epsilons):
    ax3.plot(
        rhos,
        [2 * x for x in data[(-1, eps)]["all"]["sd1"]],
        markers[ii] + "-",
        color=colors[ii],
        ms=4,
        lw=0.5,
        label=rf"$\epsilon={eps}$",
    )

for ii, eps in enumerate(epsilons):
    ax4.plot(
        rhos,
        [2 * x for x in data[(-10, eps)]["all"]["sd1"]],
        markers[ii] + "-",
        color=colors[ii],
        ms=4,
        lw=0.5,
        label=rf"$\epsilon={eps}$",
    )

ax3.set_xscale("log")
ax3.set_xlabel(r"$\rho$")
ax4.set_xscale("log")
ax4.set_xlabel(r"$\rho$")
ax3.set_ylabel(r"$\sigma_d^1$")
#ax3.legend()
#ax4.legend()

ax1.set_ylim([.005, 0.5])
ax2.set_ylim(ax2.get_ylim())
ax3.set_ylim(-1.2, 1.2)
ax4.set_ylim(ax3.get_ylim())

fig.tight_layout()
#fig.subplots_adjust(hspace=0.5, top=0.98)

#fig.text(0.05, 0.97, "A", fontsize=8, ha="center", va="center")
#fig.text(0.05, 0.45, "B", fontsize=8, ha="center", va="center")

plt.savefig(f"fig2_n0_{n0}.pdf")
