# Using the folded SFS in MSL, we fit a demographic history to the synonymous
# data, which has an ancient expansion, and recent exponential growth.
# We then fit DFEs to both missense mutations and LOF mutations. For missense
# mutations, we assume h=0.2, h=0.5, and for LOF, we assume h=0, 0.2, 0.5.


## import modules and data

import moments
import pickle
import numpy as np
import scipy.stats

data = pickle.load(open("parsed_data/MSL.frequency_spectra.folded.bp", "rb"))

fs_syn = data["synonymous"]["all"].fold()
fs_mis = data["missense"]["all"].fold()
fs_lof = data["loss_of_function"]["all"].fold()

## total mutation rates from Karczewski et al (GNOMAD mutation model)
# note that these values are u * L, where L is the total length of coding
# regions across autosomes
u_syn = 0.1442
u_mis = 0.3426
u_lof = 0.0256


## Demography inference


def model_func(params, ns):
    nuA, nuF, TA, TF = params
    fs = moments.Demographics1D.snm(ns)
    fs.integrate([nuA], TA)
    nu_func = lambda t: [nuA * np.exp(np.log(nuF / nuA) * t / TF)]
    fs.integrate(nu_func, TF)
    return fs


p_guess = [2.0, 10.0, 0.3, 0.01]
lower_bound = [1e-3, 1e-3, 1e-3, 1e-3]
upper_bound = [10, 100, 1, 1]
p_guess = moments.Misc.perturb_params(
    p_guess, lower_bound=lower_bound, upper_bound=upper_bound
)


demo_params = moments.Inference.optimize_log_fmin(
    p_guess, fs_syn, model_func, lower_bound=lower_bound, upper_bound=upper_bound
)


model = model_func(demo_params, fs_syn.sample_sizes)
opt_syn_theta = moments.Inference.optimal_sfs_scaling(model, fs_syn)
Ne = opt_syn_theta / u_syn / 4

print("optimal demog. parameters:", [float(f"{o:.4f}") for o in demo_params])
print("LL:", f"{moments.Inference.ll_multinom(model, fs_syn):0.1f}")
print("inferred Ne:", f"{Ne:.2f}")


## Cache selection spectra for h=0, h=0.2, h=0.5


def selection_spectrum(gamma, h=0.5):
    rerun = True
    ns_sim = 100
    while rerun:
        ns_sim = 2 * ns_sim
        fs = moments.LinearSystem_1D.steady_state_1D(ns_sim, gamma=gamma, h=h)
        fs = moments.Spectrum(fs)
        fs.integrate([demo_params[0]], demo_params[2], gamma=gamma, h=h)
        nu_func = lambda t: [
            demo_params[0]
            * np.exp(np.log(demo_params[1] / demo_params[0]) * t / demo_params[3])
        ]
        fs.integrate(nu_func, demo_params[3], gamma=gamma, h=h)
        if abs(np.max(fs)) > 10 or np.any(np.isnan(fs)):
            # large gamma-values can require large sample sizes for stability
            rerun = True
        else:
            rerun = False
    fs = fs.project(fs_syn.sample_sizes)
    return fs


def build_cache(gammas, h):
    spectrum_cache = {}
    spectrum_cache[0] = selection_spectrum(0)
    for gamma in gammas:
        spectrum_cache[gamma] = selection_spectrum(-gamma, h=h)
    return spectrum_cache


gammas = np.logspace(-4, 3, 101)
try:
    caches = pickle.load(open("data/supp/sfs_inference/MSL_cache.bp", "rb"))
except IOError:
    caches = {}
    for h in [0, 0.2, 0.5]:
        caches[h] = build_cache(gammas, h=h)
    pickle.dump(caches, open("data/supp/sfs_inference/MSL_cache.bp", "wb+"))


## Fit DFEs


def raw_sel_coeff_bins(shape, scale, Ne, bins=[1e-5, 1e-4, 1e-3, 1e-2]):
    ps = [scipy.stats.gamma.cdf(2 * Ne * bins[0], shape, scale=scale)]
    for ii, s in enumerate(bins[1:]):
        ps.append(
            scipy.stats.gamma.cdf(2 * Ne * s, shape, scale=scale)
            - scipy.stats.gamma.cdf(2 * Ne * bins[ii], shape, scale=scale)
        )
    ps.append(1 - scipy.stats.gamma.cdf(2 * Ne * bins[-1], shape, scale=scale))
    return ps


theta_mis = opt_syn_theta * u_mis / u_syn
theta_lof = opt_syn_theta * u_lof / u_syn

dxs = (gammas - np.concatenate(([gammas[0]], gammas))[:-1]) / 2 + (
    np.concatenate((gammas, [gammas[-1]]))[1:] - gammas
) / 2


def dfe_func(params, ns, theta=1, h=0.5):
    alpha, beta = params
    fs = caches[h][0] * scipy.stats.gamma.cdf(gammas[0], alpha, scale=beta)
    weights = scipy.stats.gamma.pdf(gammas, alpha, scale=beta)
    for gamma, dx, w in zip(gammas, dxs, weights):
        fs += caches[h][gamma] * dx * w
    fs = theta * fs
    return fs


print("Missense DFEs:")
print()

fout = open("../tables/msl_dfes.txt", "w+")
fout.write(
    "Class;$h$;shape;scale;LL;$[0,10^{-5})$;$[10^{-5},10^{-4})$;$[10^{-4},10^{-3})$;$[10^{-3},10^{-2})$;$[10^{-2},\infty)$\n"
)

for h in [0.0, 0.2, 0.5]:
    print()

    def model_func_missense(params, ns):
        return dfe_func(params, ns, theta=theta_mis, h=h)

    p_guess = [0.2, 1000]
    lower_bound = [1e-4, 1e-1]
    upper_bound = [1e1, 1e8]
    p_guess = moments.Misc.perturb_params(
        p_guess, lower_bound=lower_bound, upper_bound=upper_bound
    )

    opt_params_mis = moments.Inference.optimize_log_fmin(
        p_guess,
        fs_mis,
        model_func_missense,
        lower_bound=lower_bound,
        upper_bound=upper_bound,
        multinom=False,
    )

    model_mis = model_func_missense(opt_params_mis, fs_mis.sample_sizes)
    shape, scale = opt_params_mis
    print(f"optimal parameters, h={h}:")
    print("shape:", f"{shape:.3f}")
    print("scale:", f"{scale:.1f}")
    LL = moments.Inference.ll(model_mis, fs_mis)
    print("LL:", f"{LL:0.1f}")
    ps = raw_sel_coeff_bins(opt_params_mis[0], opt_params_mis[1], Ne)
    print("proportions in bins:", [float(f"{p:0.3f}") for p in ps])
    if h == 0.0:
        c = "Missense"
    else:
        c = ""
    fout.write(
        f"{c};{h};{shape:0.3f};{scale:0.0f};{LL:0.1f};"
        + ";".join([f"{p:0.3f}" for p in ps])
        + "\n"
    )

print("Loss-of-function DFEs:")
print()

for h in [0.0, 0.2, 0.5]:
    print()

    def model_func_lof(params, ns):
        return dfe_func(params, ns, theta=theta_lof, h=h)

    p_guess = [0.2, 1000]
    lower_bound = [1e-4, 1e-1]
    upper_bound = [1e1, 1e8]
    p_guess = moments.Misc.perturb_params(
        p_guess, lower_bound=lower_bound, upper_bound=upper_bound
    )

    opt_params_lof = moments.Inference.optimize_log_fmin(
        p_guess,
        fs_lof,
        model_func_lof,
        lower_bound=lower_bound,
        upper_bound=upper_bound,
        multinom=False,
    )

    model_lof = model_func_lof(opt_params_lof, fs_lof.sample_sizes)
    shape, scale = opt_params_lof
    print(f"optimal parameters, h={h}:")
    print("shape:", f"{shape:.3f}")
    print("scale:", f"{scale:.1f}")
    LL = moments.Inference.ll(model_lof, fs_lof)
    print("LL:", f"{LL:0.1f}")
    ps = raw_sel_coeff_bins(shape, scale, Ne)
    print("proportions in bins:", [float(f"{p:0.3f}") for p in ps])
    if h == 0.0:
        c = "LOF"
    else:
        c = ""
    fout.write(
        f"{c};{h};{shape:0.3f};{scale:0.0f};{LL:0.1f};"
        + ";".join([f"{p:0.3f}" for p in ps])
        + "\n"
    )

fout.close()

## plot demography and fits

import matplotlib.pylab as plt

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

fig = plt.figure(983701, figsize=(6.5, 6.5))
fig.clf()

ax1 = plt.subplot2grid((4, 2), (0, 0), rowspan=2)

T1 = (demo_params[2] + demo_params[3]) * 2 * Ne * 29 / 1000
T2 = demo_params[3] * 2 * Ne * 29 / 1000
ax1.plot((500, T1), (Ne, Ne), '-', color=colors[0], lw=3)
ax1.plot((T1, T1), (Ne, demo_params[0] * Ne), '-', color=colors[0], lw=3)
ax1.plot((T1, T2), (demo_params[0] * Ne, demo_params[0] * Ne), '-', color=colors[0], lw=3)
xx = np.linspace(T2, 0, 101)
yy = Ne * demo_params[0] * np.exp(np.log(demo_params[1] / demo_params[0]) * np.linspace(0, demo_params[3], 101) / demo_params[3])
ax1.plot(xx, yy, '-', color=colors[0], lw=3)

ax1.set_xlim([500, 0])
ax1.set_xlabel("Time ago (ky)")
ax1.set_ylabel("Effective population size")
ax1.set_ylim(bottom=0)

# plot syn fit and residuals

ax2 = plt.subplot2grid((4, 2), (0, 1))
ax3 = plt.subplot2grid((4, 2), (1, 1))

ax2.plot(fs_syn, "ro", label="Syn. data", ms=3, lw=0.5, fillstyle="none", alpha=0.7)
ax2.plot(opt_syn_theta * model.fold(), "k.--", label="Model", ms=2, lw=0.5, alpha=0.7)

ax2.set_yscale("log")
ax2.set_xlabel("Allele count")
ax2.set_ylabel("Num. mutations")
ax2.legend()
ax2.set_xlim(0, fs_syn.sample_sizes[0] / 2 + 1)
ax2.set_title("Demography fit to synonymous data")

resid = moments.Inference.Anscombe_Poisson_residual(opt_syn_theta * model, fs_syn)

ax3.plot(resid, ".-", lw=0.5, color=colors[3], label="Anscombe residual")
ax3.plot((0, fs_syn.sample_sizes), (0, 0), "k--", lw=0.5, label=None)
ax3.set_xlim(ax2.get_xlim())

ax3.set_xlabel("Allele count")
ax3.set_ylabel("Residual")

# plot mis fit and residuals

ax4 = plt.subplot2grid((4, 2), (2, 0))
ax5 = plt.subplot2grid((4, 2), (3, 0))

ax4.plot(fs_mis, "ro", label="Mis. data", ms=3, lw=0.5, fillstyle="none", alpha=0.7)
ax4.plot(model_mis.fold(), "k.--", label="Model", ms=2, lw=0.5, alpha=0.7)

ax4.set_yscale("log")
ax4.set_xlabel("Allele count")
ax4.set_ylabel("Num. mutations")
ax4.legend()
ax4.set_xlim(ax2.get_xlim())
ax4.set_title("DFE fit to missense data")

resid = moments.Inference.Anscombe_Poisson_residual(model_mis, fs_mis)

ax5.plot(resid, ".-", lw=0.5, color=colors[3], label="Anscombe residual")
ax5.plot((0, fs_syn.sample_sizes), (0, 0), "k--", lw=0.5, label=None)
ax5.set_xlim(ax2.get_xlim())

ax5.set_xlabel("Allele count")
ax5.set_ylabel("Residual")

# plot lof fit and residuals

ax6 = plt.subplot2grid((4, 2), (2, 1))
ax7 = plt.subplot2grid((4, 2), (3, 1))

ax6.plot(fs_lof, "ro", label="LOF data", ms=3, lw=0.5, fillstyle="none", alpha=0.7)
ax6.plot(model_lof.fold(), "k.--", label="Model", ms=2, lw=0.5, alpha=0.7)

ax6.set_yscale("log")
ax6.set_xlabel("Allele count")
ax6.set_ylabel("Num. mutations")
ax6.legend()
ax6.set_xlim(0, fs_syn.sample_sizes[0] / 2 + 1)
ax6.set_title("DFE fit to loss-of-function data")

resid = moments.Inference.Anscombe_Poisson_residual(model_lof, fs_lof)

ax7.plot(resid, ".-", lw=0.5, color=colors[3], label="Anscombe residual")
ax7.plot((0, fs_syn.sample_sizes), (0, 0), "k--", lw=0.5, label=None)
ax7.set_xlim(ax2.get_xlim())

ax7.set_xlabel("Allele count")
ax7.set_ylabel("Residual")

fig.tight_layout()
plt.savefig("../figures/msl_demography_dfes.pdf")


