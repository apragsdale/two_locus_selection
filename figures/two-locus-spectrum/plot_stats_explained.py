# Unconditioned vs allele count-conditioned LD.  Four panels, first two showing
# Hudson slices in the with rho=1 and rho=10 for neutrality and with gamma=-5
# at both loci. Second two panels show \sigma_d^2 and \sigma_d^1 for nuetrality
# and using all SNPs with gamma=-5, as well as conditioning on doubletons,
# tripletons, and few other frequencies.

import moments, moments.TwoLocus
import matplotlib.pylab as plt
import numpy as np
import pickle, gzip
import copy

# set font sizes
import matplotlib

plt.rcParams["legend.title_fontsize"] = "xx-small"
matplotlib.rc("xtick", labelsize=6)
matplotlib.rc("ytick", labelsize=6)
matplotlib.rc("axes", labelsize=7)
matplotlib.rc("axes", titlesize=7)
matplotlib.rc("legend", fontsize=6)

plt.rc("text", usetex=True)

from bokeh.palettes import Colorblind

colors = Colorblind[8]

import mpl_toolkits.mplot3d as mplot3d

# adapted from dadi's Plotting methods
def plot_3d(ax, F, vmin=None, vmax=None, max_sum=None, max_index=None):
    """
    vmin - minimum value to plot
    vmax - maximum value to plot
    max_sum - slice spectrum on diagonal plane, plotting points closer to the origin than the plane
    max_index - only show entries this distance from one of the axes planes
    """
    if vmin == None:
        vmin = np.min(F[F > 0])
    if vmax == None:
        vmax = np.max(F[F > 0])

    F2 = copy.deepcopy(F)
    F2 = np.swapaxes(F2, 0, 2)

    toplot = np.logical_not(F2.mask)
    toplot = np.logical_and(toplot, F2.data >= vmin)
    if max_sum != None:
        iis = np.where(toplot)[0]
        jjs = np.where(toplot)[1]
        kks = np.where(toplot)[2]
        for ll in range(len(iis)):
            ii = iis[ll]
            jj = jjs[ll]
            kk = kks[ll]
            if ii + jj + kk > max_sum:
                toplot[ii, jj, kk] = False
    if max_index != None:
        iis = np.where(toplot)[0]
        jjs = np.where(toplot)[1]
        kks = np.where(toplot)[2]
        for ll in range(len(iis)):
            ii = iis[ll]
            jj = jjs[ll]
            kk = kks[ll]
            if ii > max_index and jj > max_index and kk > max_index:
                toplot[ii, jj, kk] = False

    normalized = (np.log10(F2) - np.log10(vmin)) / (np.log10(vmax) - np.log10(vmin))
    normalized = np.minimum(normalized, 1)
    # scrunch by a factor
    # XXX: this is really hacky
    factor = 0.1
    normalized = (1 - 2 * factor) * normalized + 0.2
    colors = plt.cm.magma_r(normalized)

    # We draw by calculating which faces are visible and including each as a
    # polygon.
    polys, polycolors = [], []
    for ii in range(F.shape[0]):
        for jj in range(F.shape[1]):
            for kk in range(F.shape[2]):
                if not toplot[ii, jj, kk]:
                    continue
                if kk < F.shape[2] - 1 and toplot[ii, jj, kk + 1]:
                    pass
                else:
                    polys.append(
                        [
                            [ii - 0.5, jj + 0.5, kk + 0.5],
                            [ii + 0.5, jj + 0.5, kk + 0.5],
                            [ii + 0.5, jj - 0.5, kk + 0.5],
                            [ii - 0.5, jj - 0.5, kk + 0.5],
                        ]
                    )
                    polycolors.append(colors[ii, jj, kk])
                if kk > 0 and toplot[ii, jj, kk - 1]:
                    pass
                else:
                    polys.append(
                        [
                            [ii - 0.5, jj + 0.5, kk - 0.5],
                            [ii + 0.5, jj + 0.5, kk - 0.5],
                            [ii + 0.5, jj - 0.5, kk - 0.5],
                            [ii - 0.5, jj - 0.5, kk - 0.5],
                        ]
                    )
                    polycolors.append(colors[ii, jj, kk])
                if jj < F.shape[1] - 1 and toplot[ii, jj + 1, kk]:
                    pass
                else:
                    polys.append(
                        [
                            [ii - 0.5, jj + 0.5, kk + 0.5],
                            [ii + 0.5, jj + 0.5, kk + 0.5],
                            [ii + 0.5, jj + 0.5, kk - 0.5],
                            [ii - 0.5, jj + 0.5, kk - 0.5],
                        ]
                    )
                    polycolors.append(colors[ii, jj, kk])
                if jj > 0 and toplot[ii, jj - 1, kk]:
                    pass
                else:
                    polys.append(
                        [
                            [ii - 0.5, jj - 0.5, kk + 0.5],
                            [ii + 0.5, jj - 0.5, kk + 0.5],
                            [ii + 0.5, jj - 0.5, kk - 0.5],
                            [ii - 0.5, jj - 0.5, kk - 0.5],
                        ]
                    )
                    polycolors.append(colors[ii, jj, kk])
                if ii < F.shape[0] - 1 and toplot[ii + 1, jj, kk]:
                    pass
                else:
                    polys.append(
                        [
                            [ii + 0.5, jj - 0.5, kk + 0.5],
                            [ii + 0.5, jj + 0.5, kk + 0.5],
                            [ii + 0.5, jj + 0.5, kk - 0.5],
                            [ii + 0.5, jj - 0.5, kk - 0.5],
                        ]
                    )
                    polycolors.append(colors[ii, jj, kk])
                if ii > 0 and toplot[ii - 1, jj, kk]:
                    pass
                else:
                    polys.append(
                        [
                            [ii - 0.5, jj - 0.5, kk + 0.5],
                            [ii - 0.5, jj + 0.5, kk + 0.5],
                            [ii - 0.5, jj + 0.5, kk - 0.5],
                            [ii - 0.5, jj - 0.5, kk - 0.5],
                        ]
                    )
                    polycolors.append(colors[ii, jj, kk])

    polycoll = mplot3d.art3d.Poly3DCollection(
        polys, facecolor=polycolors, edgecolor="k", linewidths=0.5
    )
    ax.add_collection(polycoll)

    # Set the limits
    ax.set_xlim3d(-0.5, F.shape[0] - 0.5)
    ax.set_ylim3d(-0.5, F.shape[1] - 0.5)
    ax.set_zlim3d(-0.5, F.shape[2] - 0.5)

    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])
    ax.view_init(elev=20.0, azim=-25)
    ax.set_axis_off()
    ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))


def get_r2(F):
    r2 = 0
    n = F.sample_size
    for i in range(n):
        for j in range(n):
            for k in range(n - i - j):
                nA = i + j
                if nA == n or nA == 0:
                    continue
                fA = nA / n
                nB = i + k
                if nB == n or nB == 0:
                    continue
                fB = nB / n
                fAB = i / n
                fAb = j / n
                faB = k / n
                fab = (n - i - j - k) / n
                r2 += (
                    F[i, j, k] * (fAB - fA * fB) ** 2 / (fA * (1 - fA) * fB * (1 - fB))
                )
        return r2


n = 40
gamma = -2
rho = 1

L = 100e6
theta = 4 * 1e4 * 1.25e-8 * L

fs = theta * moments.Spectrum(moments.LinearSystem_1D.steady_state_1D(n, gamma=gamma))
fs = fs.project([10])

try:
    F = pickle.load(gzip.open("tls.fs", "rb"))
except IOError:
    F = moments.TwoLocus.Demographics.equilibrium(n, gamma=gamma, rho=rho)
    F = F.project(10)
    pickle.dump(F, gzip.open("tls.fs", "wb+"))

try:
    rhos, sd1, sd2 = pickle.load(open("sigma_data.bp", "rb"))
except IOError:
    rhos = np.logspace(-2, 1.2, 21)
    for rho in rhos:
        F_sel = moments.TwoLocus.Demographics.equilibrium(30, gamma=gamma, rho=rho)
        sd1.append(F_sel.D() / F_sel.pi2())
        sd2.append(F_sel.D2() / F_sel.pi2())
    pickle.dump([rhos, sd1, sd2], open("sigma_data.bp", "wb+"))

fs = fs.project([10])
F = F.project(10)

fig = plt.figure(1, figsize=(6.5, 3.5))
fig.clf()

grids = (2, 10)

ax1 = plt.subplot2grid(grids, (0, 0), colspan=4)

fs_neu = moments.Demographics1D.snm([10]) * theta
ax1.plot(fs_neu, "o--", color=colors[0], ms=3, lw=0.5, label="Neutral")
ax1.plot(fs, "o--", color=colors[1], ms=3, lw=0.5, label="Selected")
ax1.set_xlabel("Derived allele count")
ax1.set_ylabel(r"Mutation count")
ax1.set_xlim([0, 10])
ax1.set_ylim(bottom=0)
ax1.spines["top"].set_visible(False)
ax1.spines["right"].set_visible(False)

# ax1.set_yscale("log")

# fig.text(
#    0.30,
#    0.95,
ax1.set_title("Single-locus site frequency spectrum ($2Ns = -2$)", fontsize=7)
#    va="center",
#    ha="center",
# )

# inset panel
ax1.text(3, 50000, "\\underline{Summaries (selected):}", fontsize=6, ha="left", va="center")
ax1.text(
    3, 44000, rf"$\pi$ per bp = ${fs.pi() / L:.5f}$", fontsize=6, ha="left", va="center"
)
ax1.text(
    3, 38000, rf"Tajima's $D = {fs.Tajima_D():.3f}$", fontsize=6, ha="left", va="center"
)

ax1.legend(loc="center right")

ax2 = plt.subplot2grid(grids, (0, 4), colspan=4, projection="3d")

plot_3d(ax2, F)

ax2.text(0, -2, 0, "0", fontsize=6)
ax2.text(0, -2.5, 10, "10", fontsize=6)
ax2.text(0, -6, 5, "$n_{AB}$")

ax2.text(6, -5, 0, "$n_{Ab}$")
ax2.text(1, -2.5, -1, "0", fontsize=6)
ax2.text(10, -2.5, -1, "10", fontsize=6)

ax2.text(0, 12, 0, "$n_{aB}$")

fig.text(
    0.70,
    0.95,
    "Two-locus haplotype sampling distribution",
    fontsize=7,
    va="center",
    ha="center",
)
fig.text(
    0.70,
    0.91,
    r"($2Ns_A = -2, 2Ns_B=-2, \epsilon=0$)",
    fontsize=7,
    va="center",
    ha="center",
)

fig.text(0.75, 0.80, "\\underline{Summaries:}", fontsize=6, ha="left", va="center")
fig.text(0.75, 0.75, rf"$r^2 = {get_r2(F):.3f}$", fontsize=6, ha="left", va="center")
fig.text(
    0.75,
    0.70,
    rf"$\sigma_d^2 = {F.D2() / F.pi2():.3f}$",
    fontsize=6,
    ha="left",
    va="center",
)
fig.text(
    0.75,
    0.65,
    rf"$\sigma_d^1 = {F.D() / F.pi2():.3f}$",
    fontsize=6,
    ha="left",
    va="center",
)

ax3 = plt.subplot2grid(grids, (1, 0), colspan=4)
ax3.set_title("$n_{AB}$ proportions with $n_A=5, n_B=5$", fontsize=7)
ax3.set_ylabel("Proportion")
ax3.set_xlabel("$n_{AB}$")
c, p = moments.TwoLocus.Util.pAB(F, 5, 5)
F_neu = moments.TwoLocus.Demographics.equilibrium(30, rho=1).project(10)
c, p_neu = moments.TwoLocus.Util.pAB(F_neu, 5, 5)
ax3.bar(c-0.2, p_neu / p_neu.sum(), width=0.35, label="$2Ns=0$", color=colors[0])
ax3.bar(c+0.2, p / p.sum(), width=0.35, label="$2Ns=-2$", color=colors[1])
ax3.legend()

ax4 = plt.subplot2grid(grids, (1, 4), colspan=3)
ax4.set_title("Decay of $\sigma_d^2$ with recomb. distance", fontsize=7)
ax4.set_ylabel(r"$\sigma_d^2$")
ax4.set_xlabel(r"$\rho=4Nr$")
rhos = np.logspace(-2, 1.2, 21)
ok = (5 + rhos / 2) / (11 + 13 * rhos / 2 + rhos ** 2 / 2)
ax4.plot(rhos, ok, "-", color=colors[0], lw=1, label=r"$2Ns=0$ (Ohta \& Kimura)")
ax4.plot(rhos, sd2, "-", color=colors[1], ms=3, lw=1, label=r"$2Ns=-2$")
ax4.set_xscale("log")
ax4.set_yscale("log")
ax4.set_ylim(0.01, 1)
ax4.legend(fontsize=5)

ax5 = plt.subplot2grid(grids, (1, 7), colspan=3)
ax5.set_title("$\sigma_d^1$ against recomb. distance", fontsize=7)
ax5.set_xlabel(r"$\rho=4Nr$")
ax5.set_ylabel(r"$\sigma_d^1$")
ok = 0 * rhos
ax5.plot(rhos, ok, "-", color=colors[0], lw=1, label=r"$2Ns=0$")
ax5.plot(rhos, sd1, "-", color=colors[1], ms=3, lw=1, label=r"$2Ns=-2$")
ax5.set_xscale("log")
ax5.set_ylim(-0.3, 0.1)
ax5.legend(fontsize=5)

fig.tight_layout()
fig.text(0.02, 0.95, "A", fontsize=8, ha="center", va="center")
fig.text(0.52, 0.95, "B", fontsize=8, ha="center", va="center")
fig.text(0.04, 0.47, "C", fontsize=8, ha="center", va="center")
fig.text(0.43, 0.47, "D", fontsize=8, ha="center", va="center")
fig.text(0.73, 0.47, "E", fontsize=8, ha="center", va="center")
plt.savefig("stats_summaries.pdf")
