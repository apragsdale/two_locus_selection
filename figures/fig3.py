# Figure 3: Dominance selection model.
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

n0 = int(sys.argv[1])  # larger sample size for improved accuracy
n = 30  # the projection size

hs = [0.0, 0.1, 0.5, 1.0]
gammas = [-1, -5]
conditions = ["all", 2, 3, 4, 10]

data_fname = f"fig3_data_ns_{n0}_{n}.bp"
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


fig = plt.figure(3, figsize=(6.5, 3.5))
fig.clf()

rhos = data["rhos"]
ok = data["Ohta_Kimura"]

markers = ["x", "+", "1", "2"]
colors = Colorblind[8]

ax1 = plt.subplot(2, 2, 1)
ax2 = plt.subplot(2, 2, 2)

ax1.plot(rhos, ok, "k--", lw=1, label=None)
ax2.plot(rhos, ok, "k--", lw=1, label=None)

for ii, h in enumerate(hs):
    ax1.plot(
        rhos,
        data[(-1, h)]["all"]["sd2"],
        markers[ii] + "-",
        color=colors[ii],
        ms=4,
        lw=0.5,
        label=rf"$h={h}$",
    )

for ii, h in enumerate(hs):
    ax2.plot(
        rhos,
        data[(-5, h)]["all"]["sd2"],
        markers[ii] + "-",
        color=colors[ii],
        ms=4,
        lw=0.5,
        label=rf"$h={h}$",
    )

ax3 = plt.subplot(2, 2, 3)
ax4 = plt.subplot(2, 2, 4)
ax3.plot(rhos, 0 * ok, "k--", lw=1, label=None)
ax4.plot(rhos, 0 * ok, "k--", lw=1, label=None)

for ii, h in enumerate(hs):
    ax3.plot(
        rhos,
        [2 * x for x in data[(-1, h)]["all"]["sd1"]],
        markers[ii] + "-",
        color=colors[ii],
        ms=4,
        lw=0.5,
        label=rf"$h={h}$",
    )

for ii, h in enumerate(hs):
    ax4.plot(
        rhos,
        [2 * x for x in data[(-5, h)]["all"]["sd1"]],
        markers[ii] + "-",
        color=colors[ii],
        ms=4,
        lw=0.5,
        label=rf"$h={h}$",
    )

ax1.set_xscale("log")
ax1.set_yscale("log")
ax2.set_xscale("log")
ax2.set_yscale("log")
ax3.set_xscale("log")
ax4.set_xscale("log")
ax1.set_ylabel(r"$\sigma_d^2$")
ax3.set_ylabel(r"$\sigma_d^1$")
ax3.set_xlabel(r"$\rho$")
ax4.set_xlabel(r"$\rho$")
ax1.legend()

ax1.set_title(r"$\gamma=-1$")
ax2.set_title(r"$\gamma=-5$")

fig.tight_layout()
#fig.text(0.05, 0.97, "A", fontsize=8, ha="center", va="center")
#fig.text(0.05, 0.47, "B", fontsize=8, ha="center", va="center")

plt.savefig(f"fig3_n0_{n0}.pdf")
