import gzip, os, pickle
import moments.TwoLocus
import numpy as np
import matplotlib.pylab as plt

N = 1000
n = 50

hs = np.linspace(0.05, 0.5, 10)
rs = [0.001, 0.0001]
shs = [0.01, 0.001, 0.0004, 0.0001]

#files = os.listdir("outputs")
#roze_simulations = []
#for f in files:
#    if "_sh_" in f:
#        roze_simulations.append(f)
#
#vals = {}
#for f in roze_simulations:
#    data = pickle.load(gzip.open("outputs/" + f, "rb"))
#    f = f.split(".bp")[0]
#    r = float(f.split("_r_")[1].split("_")[0])
#    sh = float(f.split("_sh_")[1].split("_")[0])
#    h = float(f.split("_h_")[1].split("_")[0])
#    F = data["simulation"]
#    vals.setdefault(r, {})
#    vals[r].setdefault(sh, {})
#    vals[r][sh][h] = F.D() / F.pi2()

vals = pickle.load(open("roze_sd1s.bp", "rb"))

#exp = {}
#rhos = [4 * N * r for r in rs]
#spectra = os.listdir("data")
#for f in spectra:
#    if "sh" in f:
#        F = pickle.load(gzip.open("data/" + f, "rb"))
#        f = f.split(".bp")[0]
#        rho = float(f.split("_rho_")[1].split("_")[0])
#        sh = float(f.split("_sh_")[1].split("_")[0])
#        h = float(f.split("_h_")[1].split("_")[0])
#        exp.setdefault(rho, {})
#        exp[rho].setdefault(sh, {})
#        exp[rho][sh][h] = F.D() / F.pi2()

exp = pickle.load(open("roze_comp_data.bp", "rb"))

fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(6, 3))
r = 0.0001
for sh in shs[::-1]:
    v = [vals[r][sh][h] for h in hs]
    ax1.plot(hs, v, "o--", ms=3, lw=0.5, label=f"{sh} (sim.)")

ax1.set_prop_cycle(None)
for sh in shs[::-1]:
    v = [exp[4 * N * r][sh][h] for h in sorted(exp[4 * N * r][sh].keys())]
    if sh == 0.01:
        continue
    elif sh == 0.001:
        ax1.plot(hs[1:], v[1:], "x-", ms=5, lw=1, label="moments")
    else:
        ax1.plot(hs, v, "x-", ms=5, lw=1, label="moments")

ax1.set_title(f"$r={r}$")
ax1.legend(title="$sh$", ncol=2, fontsize=6, loc="lower right")
ax1.set_xlabel("$h$")
ax1.set_xlim(0, 0.55)
ax1.set_ylim(-0.55, 0.15)

r = 0.001
for sh in shs[::-1]:
    v = [vals[r][sh][h] for h in hs]
    ax2.plot(hs, v, "o--", ms=3, lw=0.5, label=f"{sh} (sim.)")

ax2.set_prop_cycle(None)
for sh in shs[::-1]:
    v = [exp[4 * N * r][sh][h] for h in sorted(exp[4 * N * r][sh].keys())]
    if sh == 0.01:
        continue
    elif sh == 0.001:
        ax2.plot(hs[1:], v[1:], "x-", ms=5, lw=1, label="moments")
    else:
        ax2.plot(hs, v, "x-", ms=5, lw=1, label="moments")

ax2.set_title(f"$r={r}$")
ax2.legend(title="$sh$", ncol=2, fontsize=6, loc="lower right")
ax2.set_xlabel("$h$")
ax2.set_xlim(0, 0.55)
ax2.set_ylim(-0.125, 0.075)
ax1.set_ylabel("$\sigma_d^1$")
fig.tight_layout()
plt.savefig("Roze.pdf")
