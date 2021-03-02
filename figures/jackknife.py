import moments.TwoLocus
import numpy as np
import matplotlib.pylab as plt

import matplotlib
plt.rcParams["legend.title_fontsize"] = "xx-small"
matplotlib.rc("xtick", labelsize=6)
matplotlib.rc("ytick", labelsize=6)
matplotlib.rc("axes", labelsize=7)
matplotlib.rc("axes", titlesize=7)
matplotlib.rc("legend", fontsize=6)

from matplotlib.ticker import ScalarFormatter
ns = [10, 15, 20, 25, 30, 35, 40, 50, 60, 70]
t_jk = np.array(
    [1.01, 1.48, 3.03, 7.65, 20.58, 47.07, 102.44, 412.94, 1192.44, 3526.23]
)
mem_jk = np.array(
    [131828, 135580, 141608, 157068, 187996, 253480, 345692, 700044, 1601448, 3543684]
)
t_no_jk = np.array(
    [0.87, 0.93, 1.02, 1.43, 2.21, 4.24, 8.31, 26.66, 77.12, 280.66]
)
mem_no_jk = np.array(
    [130864, 133832, 140784, 156116, 185932, 247812, 341120, 692348, 1598728, 3541400]
)



rhos_ok = np.logspace(-1, 2, 50)
ohta_kimura = (5 + rhos_ok / 2) / (11 + 13 * rhos_ok / 2 + rhos_ok ** 2 / 2)
rhos = np.logspace(-1, 2, 11)
ok_compare = (5 + rhos / 2) / (11 + 13 * rhos / 2 + rhos ** 2 / 2)
# precomputed using `F = moments.TwoLocus.Demographics.equilibrium(n, rho=rho)`
# and then `F.D2() / F.pi2()`
ld_curve_moments = {
    20: [0.4332, 0.4138, 0.3799, 0.3264, 0.2547, 0.1774, 0.1108, 0.0634, 0.0339, 0.0172, 0.0045],
    30: [0.4332, 0.4139, 0.3801, 0.3269, 0.2556, 0.1786, 0.1121, 0.0646, 0.035, 0.0248, 0.0074],
    50: [0.4333, 0.414, 0.3803, 0.3272, 0.2562, 0.1794, 0.1128, 0.0652, 0.0356, 0.0186, 0.2883],
    80: [0.4333, 0.414, 0.3803, 0.3273, 0.2565, 0.1797, 0.1131, 0.0655, 0.0357, -0.0117, -0.6302],
}

lw = 1
ms = 3

fig = plt.figure(figsize=(6.5, 5))
ax1 = plt.subplot(2, 2, 1)
ax1.plot(rhos_ok, ohta_kimura, 'k--', lw=lw, label="Ohta and Kimura")
for n in sorted(ld_curve_moments.keys()):
    ax1.plot(rhos, ld_curve_moments[n], "v-", lw=lw, ms=ms, label=f"moments, n={n}")
ax1.set_ylabel(r"$\sigma_d^2$")
ax1.set_yscale("log")
ax1.set_xscale("log")
ax1.legend()
ax1.set_title("Comparison to analytic expectation")
ax1.set_xlabel(r"$\rho$")

ax2 = plt.subplot(2, 2, 2)
ax2.plot(rhos_ok, rhos_ok * 0, "k--", lw=lw, label=None)
for n in sorted(ld_curve_moments.keys()):
    ax2.plot(
        rhos[:-2],
        ((ld_curve_moments[n] - ok_compare) / ok_compare)[:-2],
        "v-",
        lw=1,
        ms=ms,
        label=f"moments, n={n}"
    )
ax2.set_xlim(ax1.get_xlim())
ax2.set_ylabel("Percent error")
ax2.set_xlabel(r"$\rho$")
ax2.set_xscale("log")
ax2.legend()
ax2.set_title("Moment closure accuracy")

ax3 = plt.subplot(2, 2, 3)
ax3.plot(ns, t_jk, lw=lw, label="With jackknife computation")
ax3.plot(ns, t_no_jk, lw=lw, label="Cached jackknife")
ax3.set_xlabel("Sample size")
ax3.set_ylabel("Time (seconds)")
ax3.set_yscale("log")
ax3.set_xscale("log")
ax3.legend()
ax3.xaxis.set_major_formatter(ScalarFormatter())
ax3.xaxis.set_minor_formatter(ScalarFormatter())
ax3.yaxis.set_major_formatter(ScalarFormatter())
ax3.set_title("Time to compute equilibrium FS")

ax4 = plt.subplot(2, 2, 4)
ax4.plot(ns, mem_jk / 1024, lw=lw, label="With jackknife computation")
ax4.plot(ns, mem_no_jk / 1024, lw=lw, label="Cached jackknife")
ax4.set_xlabel("Sample size")
ax4.set_ylabel("Mb")
ax4.set_yscale("log")
ax4.set_xscale("log")
ax4.legend()
ax4.xaxis.set_major_formatter(ScalarFormatter())
ax4.xaxis.set_minor_formatter(ScalarFormatter())
ax4.yaxis.set_major_formatter(ScalarFormatter())
ax4.set_title("Maximum memory usage")

fig.tight_layout()
plt.savefig("jackknife.pdf")
