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

gammas = -np.concatenate((np.logspace(-1, 1, 21), [15.0, 20.0]))

rho = 0.0
hs = [0.05, 0.25, 0.45]

data = {}
for h in hs:
    print("h =", h)
    data[h] = []
    for gamma in gammas:
        print("gamma =", gamma, end="")
        if -gamma < 7.5:
            n = 50
        else:
            n = 80
        sel_params = moments.TwoLocus.Util.simple_dominance(gamma, h=h)
        F = moments.TwoLocus.Demographics.equilibrium(
            n, sel_params_general=sel_params, rho=rho
        )
        data[h].append(F)
        print(", min=", F.min())
    print()

sd1s = {}
for h in hs:
    print(h)
    sd1s[h] = []
    for F in data[h]:
        sd1s[h].append(F.D() / F.pi2())
    print()

import pickle, gzip

pickle.dump(
    {"spectra": data, "sd1": sd1s, "gammas": gammas},
    gzip.open("data/dominance_vary_gamma.bp.gz", "wb+"),
)
