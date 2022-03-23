# Unconditioned vs allele count-conditioned LD.  Four panels, first two showing
# Hudson slices in the with rho=1 and rho=10 for neutrality and with gamma=-5
# at both loci. Second two panels show \sigma_d^2 and \sigma_d^1 for nuetrality
# and using all SNPs with gamma=-5, as well as conditioning on doubletons,
# tripletons, and few other frequencies.

import moments.TwoLocus
import matplotlib.pylab as plt
import numpy as np
import pickle
import time
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

gammas = [-1, -20]
h = 0.05
eps = -0.95


try:
    spectra_compens = pickle.load(gzip.open("data/compensatory.bp", "rb"))
except IOError:
    rhos = np.concatenate((np.logspace(-2, 1, 22), [15, 20, 30]))
    spectra_compens = {}
    spectra_compens["rhos"] = rhos
    n = 50
    for gamma in gammas:
        spectra_compens[(gamma, eps)] = {}
        sel_params = moments.TwoLocus.Util.additive_epistasis(gamma, eps)
        for rho in rhos:
            time1 = time.time()
            print("rho:", rho)
            F_c = moments.TwoLocus.Demographics.equilibrium(
                n, rho=rho, sel_params=sel_params
            )
            spectra_compens[(gamma, eps)][rho] = F_c
            print(" finished in", int(time.time() - time1), "seconds")
            print(" min:", F_c.min())
            print(" sd1:", F_c.D() / F_c.pi2())
            print()
    pickle.dump(spectra_compens, gzip.open("data/compensatory.bp", "wb+"))

try:
    spectra_gene = pickle.load(gzip.open("data/gene_based.bp", "rb"))
except IOError:
    rhos = np.concatenate((np.logspace(-2, 1, 22), [15, 20, 30]))
    spectra_gene = {}
    spectra_gene["rhos"] = rhos
    for gamma in gammas:
        if abs(gamma) >= 5:
            n = 60
        else:
            n = 50
        print("gamma:", gamma)
        spectra_gene[(gamma, eps)] = {}
        sel_params = moments.TwoLocus.Util.gene_based_dominance(gamma, h)
        for rho in rhos:
            time1 = time.time()
            print("rho:", rho)
            F_g = moments.TwoLocus.Demographics.equilibrium(
                n, rho=rho, sel_params_general=sel_params
            )
            spectra_gene[(gamma, h)][rho] = F_g
            print(" finished in", int(time.time() - time1), "seconds")
            print(" min:", F_g.min())
            print(" sd1:", F_g.D() / F_g.pi2())
            print()
    pickle.dump(spectra_gene, gzip.open("data/gene_based.bp", "wb+"))


rhos = spectra_compens["rhos"]
ok = (5 + rhos / 2) / (11 + 13 * rhos / 2 + rhos ** 2 / 2)

# compensatory, sd2
sd1_1 = []
sd1_20 = []
sd2_1 = []
sd2_20 = []
for rho in rhos:
    F = spectra_compens[(-20, -0.95)][rho]
    D = F.D()
    D2 = F.D2()
    pi2 = F.pi2()
    sd1_20.append(D / pi2)
    sd2_20.append(D2 / pi2)

for rho in rhos:
    F = spectra_compens[(-1, -0.95)][rho]
    D = F.D()
    D2 = F.D2()
    pi2 = F.pi2()
    sd1_1.append(D / pi2)
    sd2_1.append(D2 / pi2)

sd1_1g = []
sd1_20g = []
sd2_1g = []
sd2_20g = []
for rho in rhos:
    F = spectra_gene[(-20, 0.05)][rho]
    D = F.D()
    D2 = F.D2()
    pi2 = F.pi2()
    sd1_20g.append(D / pi2)
    sd2_20g.append(D2 / pi2)

for rho in rhos:
    F = spectra_gene[(-1, 0.05)][rho]
    D = F.D()
    D2 = F.D2()
    pi2 = F.pi2()
    sd1_1g.append(D / pi2)
    sd2_1g.append(D2 / pi2)

fig = plt.figure(937492, figsize=(6.5, 3.5))
fig.clf()

ax1 = plt.subplot(2, 2, 1)
ax2 = plt.subplot(2, 2, 2)
ax3 = plt.subplot(2, 2, 3)
ax4 = plt.subplot(2, 2, 4)

ax1.plot(rhos, ok, "k--", lw=1, label="Ohta & Kimura")
ax2.plot(rhos, ok, "k--", lw=1, label="Ohta & Kimura")
ax3.plot(rhos, 0 * ok, "k--", lw=1, label=None)
ax4.plot(rhos, 0 * ok, "k--", lw=1, label=None)

# compensatory mutations sigma_d2
ax1.plot(
    rhos, sd2_20, "-", color=colors[0], lw=1, label=r"$\gamma=-20, \epsilon=-0.95$"
)
ax1.plot(rhos, sd2_1, "-", color=colors[1], lw=1, label=rf"$\gamma=-1, \epsilon=-0.95$")

# compensatory mutations sigma_d1
ax3.plot(
    rhos, sd1_20, "-", color=colors[0], lw=1, label=r"$\gamma=-20, \epsilon=-0.95$"
)
ax3.plot(rhos, sd1_1, "-", color=colors[1], lw=1, label=rf"$\gamma=-1, \epsilon=-0.95$")

# gene based dominance sigma_d2
ax2.plot(
    rhos, sd2_20g, "-", color=colors[0], lw=1, label=rf"$\gamma=-20, h=0.05$"
)
ax2.plot(
    rhos, sd2_1g, "-", color=colors[1], lw=1, label=rf"$\gamma=-1, h=0.05$"
)

# gene-based dominance sigma_d1
ax4.plot(
    rhos, sd1_20g, "-", color=colors[0], lw=1, label=rf"$\gamma=-20, h=0.05$"
)
ax4.plot(
    rhos, sd1_1g, "-", color=colors[1], lw=1, label=rf"$\gamma=-1, h=0.05$"
)


ax1.set_xscale("log")
ax2.set_xscale("log")
ax3.set_xscale("log")
ax4.set_xscale("log")
ax1.set_yscale("log")
ax2.set_yscale("log")

ax1.set_ylabel(r"$\sigma_d^2$")
ax3.set_ylabel(r"$\sigma_d^1$")

ax1.legend(loc="lower left")
ax2.legend(loc="lower left")
ax1.set_title(r"Compensatory mutation model")
ax2.set_title(r"Gene-based dominance model")

ax3.set_xlabel(r"$\rho = 4Nr$")
ax4.set_xlabel(r"$\rho = 4Nr$")

#ax1.set_ylim([0.005, 0.6])
#ax2.set_ylim([-1.2, 1.2])
#ax3.set_ylim(ax1.get_ylim())
#ax4.set_ylim(ax2.get_ylim())

fig.tight_layout()
fig.text(0.04, 0.95, "A", fontsize=8, ha="center", va="center")
fig.text(0.53, 0.95, "B", fontsize=8, ha="center", va="center")
fig.text(0.04, 0.50, "C", fontsize=8, ha="center", va="center")
fig.text(0.53, 0.50, "D", fontsize=8, ha="center", va="center")

plt.savefig(f"gene-based-compensatory.pdf")
plt.savefig(f"gene-based-compensatory.png", dpi=300)
