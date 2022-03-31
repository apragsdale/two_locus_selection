# Unconditioned vs allele count-conditioned LD.  Four panels, first two showing
# Hudson slices in the with rho=1 and rho=10 for neutrality and with gamma=-5
# at both loci. Second two panels show \sigma_d^2 and \sigma_d^1 for nuetrality
# and using all SNPs with gamma=-5, as well as conditioning on doubletons,
# tripletons, and few other frequencies.

import moments.TwoLocus
import matplotlib.pylab as plt
import numpy as np
import pickle
import gzip

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

markers = ["o", "s", "v", "^"]
ms = 2

fig = plt.figure(1, figsize=(6.5, 2))
fig.clf()

ax1 = plt.subplot(1, 3, 1)

n0 = 40
n = 30
data = pickle.load(open(f"data/fig1_data_ns_{n0}_{n}.bp", "rb"))
gammas = [-0.1, -1, -5, -10]
rhos = data["C"]["rhos"]

ax1.plot(rhos, data["C"]["Ohta_Kimura"], "k--", lw=1, label="Ohta & Kimura")

for ii, gamma in enumerate(gammas):
    ax1.plot(
        rhos,
        data["C"][gamma]["all"]["sd2"],
        # markers[ii] + "-",
        "-",
        ms=ms,
        lw=1,
        label=rf"$\gamma={gamma}$",
        color=colors[ii],
    )

ax1.set_yscale("log")
ax1.set_xscale("log")
ax1.set_ylabel(r"$\sigma_d^2$")
ax1.set_title("Normalized variance of $D$ ($\sigma_d^2$)")
ax1.set_xlabel(r"$\rho=4Nr$")

ax2 = plt.subplot(1, 3, 2)

ax2.plot(rhos, 0 * rhos, "k--", lw=1, label=None)
for ii, gamma in enumerate(gammas):
    ax2.plot(
        rhos,
        data["C"][gamma]["all"]["sd1"],
        # markers[ii] + "-",
        "-",
        ms=ms,
        lw=1,
        label=rf"$\gamma={gamma}$",
        color=colors[ii],
    )

ax2.set_xscale("log")
ax2.set_ylabel(r"$\sigma_d^1$")
ax2.set_title("Normalized signed LD ($\sigma_d^1$)")
ax2.legend(fontsize=5, loc="lower right")
ax2.set_xlabel(r"$\rho=4Nr$")

ax3 = plt.subplot(1, 3, 3)
data = pickle.load(gzip.open("data/n-80_eps-0_h-0.5.bp.gz", "rb"))

gammas = np.concatenate((np.logspace(-1, 1, 21), [15, 20, 30, 40]))
rhos = [0, 1, 10]

try:
    fname = "data/n-0_eps-0_h-0.5.sd1s.bp"
    sd1s = pickle.load(open(fname, "rb"))
except IOError:
    sd1s = {}
    for rho in rhos:
        sd1s[rho] = []
        for gamma in gammas:
            sd1s[rho].append(data[rho][gamma].D() / data[rho][gamma].pi2())
            print(rho, gamma)
    pickle.dump(sd1s, open(fname, "wb+"))

ax3.plot(gammas, 0 * gammas, "k--", lw=1, label=None)
for ii, rho in enumerate(rhos):
    ax3.plot(gammas, sd1s[rho], color=colors[ii + 4], lw=1, label=rf"$\rho={rho}$")

ax3.set_title("Normalized signed LD ($\sigma_d^1$)")
ax3.set_ylabel(r"$\sigma_d^1$")
ax3.set_xlabel(r"$\gamma=2Ns$")
ax3.set_xscale("log")
ax3.legend(fontsize=5, loc="lower right")
ax3.set_xticks([0.1, 1, 10])
ax3.set_xticklabels([-0.1, -1, -10])

fig.tight_layout()
fig.text(0.04, 0.93, "A", fontsize=8, ha="center", va="center")
fig.text(0.36, 0.93, "B", fontsize=8, ha="center", va="center")
fig.text(0.69, 0.93, "C", fontsize=8, ha="center", va="center")

plt.savefig(f"hill_robertson.pdf")
plt.savefig(f"hill_robertson.png", dpi=300)
