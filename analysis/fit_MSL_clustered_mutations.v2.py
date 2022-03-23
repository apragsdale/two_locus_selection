# Fit a model of clustered mutations to the MSL data, using sigma_d^1 and sigma_d^2

import os, sys
import scipy.optimize
import pickle
import numpy as np
import moments.TwoLocus

data = pickle.load(open("parsed_data_v2/MSL.unphased.MNM_data.bp", "rb"))

bin_edges = data["data"]["bin_edges"]
bin_mids = np.mean(bin_edges, axis=1)[:-2]

sd1 = np.array(data["data"]["bins"]["all_sites"]["synonymous"]["sd1"])[:-2]
sd2 = np.array(data["data"]["bins"]["all_sites"]["synonymous"]["sd2"])[:-2]
v1 = np.array([data["SEs"]["bins"]["all_sites"]["synonymous"][b]["sd1"] for b in bin_edges])[:-2]
v2 = np.array([data["SEs"]["bins"]["all_sites"]["synonymous"][b]["sd2"] for b in bin_edges])[:-2]

demo_params = [2.194, 5.145, 0.579, 0.0400]
Ne = 12314.8

n = 40


def demo_func(demo_params, rho, n=30, mnm=False):
    F = moments.TwoLocus.Integration.steady_state(n, rho=rho, clustered_mutations=mnm)
    F.integrate(demo_params[0], demo_params[2], rho=rho, clustered_mutations=mnm)
    nu_func = lambda t: demo_params[0] * np.exp(
        np.log(demo_params[1] / demo_params[0]) * t / demo_params[3]
    )
    F.integrate(nu_func, demo_params[3], rho=rho, clustered_mutations=mnm)
    return F


def cache_spectra(demo_params, Ne, bin_mids, r):
    try:
        cache = pickle.load(
            open(f"data/supp/mnm_inference/MSL.cache.n.{n}.r.{r}.bp", "rb")
        )
        for k in ["mnm", "ism"]:
            for b in bin_mids:
                if b not in cache[k].keys():
                    rho = 4 * Ne * b * r
                    print("recomputing", rho)
                    if k == "mnm":
                        F = demo_func(demo_params, rho, n=n, mnm=True)
                        cache["mnm"][b] = F
                    elif k == "ism":
                        F = demo_func(demo_params, rho, n=n, mnm=False)
                        cache["ism"][b] = F
        pickle.dump(
            cache, open(f"data/supp/mnm_inference/MSL.cache.n.{n}.r.{r}.bp", "wb+")
        )
    except:
        cache = {"mnm": {}, "ism": {}}
        rhos = 4 * Ne * bin_mids * r
        for rho, b in zip(rhos, bin_mids):
            F = demo_func(demo_params, rho, n=n, mnm=False)
            cache["ism"][b] = F
            F = demo_func(demo_params, rho, n=n, mnm=True)
            cache["mnm"][b] = F
            print(f"r={r}, finished rho={rho}")
        pickle.dump(
            cache, open(f"data/supp/mnm_inference/MSL.cache.n.{n}.r.{r}.bp", "wb+")
        )
    return cache


## pseudo-infer average r, but for a fixed r value, infer decay of mnm rate over bins


def ll_normal(data, model, variances):
    # model is Esd1, Esd2
    # data is sd1, sd2
    # variances is varsd1, varsd2
    Esd1, Esd2 = model
    sd1, sd2 = data
    v1, v2 = variances
    ll_sd1 = np.sum(
        -1 / 2 * (sd1 - Esd1) ** 2 / v1
        - 1 / 2 * np.log(2 * np.pi * len(sd1))
        - 1 / 2 * np.log(v1)
    )
    #ll_sd2 = np.sum(
    #    -1 / 2 * (sd2 - Esd2) ** 2 / v2
    #    - 1 / 2 * np.log(2 * np.pi * len(sd2))
    #    - 1 / 2 * np.log(v2)
    #)
    return ll_sd1 #+ ll_sd2


def object_func_log(log_params, *args, **kwargs):
    return object_func(np.exp(log_params), *args, **kwargs)


def optimize_log_fmin(
    p0, data, variances, model_func, verbose=0, lower_bound=None, upper_bound=None
):
    args = (
        model_func,
        data,
        variances,
        verbose,
        lower_bound,
        upper_bound,
    )
    outputs = scipy.optimize.fmin(
        object_func_log, np.log(p0), args=args, full_output=True, disp=False
    )
    xopt, fopt, iter, funcalls, warnflag = outputs
    return np.exp(xopt), fopt


_counter = 0


def object_func(
    params, model_func, data, variances, verbose=0, lower_bound=None, upper_bound=None,
):
    global _counter
    _counter += 1
    _out_of_bounds_val = 1e8
    output_stream = sys.stdout
    # Check our parameter bounds
    if lower_bound is not None:
        for pval, bound in zip(params, lower_bound):
            if bound is not None and pval < bound:
                return -_out_of_bounds_val
    if upper_bound is not None:
        for pval, bound in zip(params, upper_bound):
            if bound is not None and pval > bound:
                return -_out_of_bounds_val

    model = model_func(params)

    result = ll_normal(data, model, variances)

    if (verbose > 0) and (_counter % verbose == 0):
        param_str = "array([%s])" % (", ".join(["%- 12g" % v for v in params]))
        output_stream.write(
            "%-8i, %-12g, %s%s" % (_counter, result, param_str, os.linesep)
        )
    return -result


rs = [1e-9, 2e-9, 3e-9, 4e-9, 5e-9, 6e-9, 7e-9, 1e-8, 2e-8]

fout = open("mnm_optimization.v2.sd1.txt", "w+")
fout.write("r\tA\tlam\tLL\n")

for r in rs:
    print()
    print()
    print("Starting r:", r)
    cache = cache_spectra(demo_params, Ne, bin_mids, r)
    stats = {}
    for mut_type in cache.keys():
        stats[mut_type] = {}
        for b, F in cache[mut_type].items():
            stats[mut_type][b] = {"D": F.D(), "D2": F.D2(), "pi2": F.pi2()}

    def model_func(params):
        A, lam = params
        weights = A * np.exp(-lam * bin_mids)
        Esd1 = []
        Esd2 = []
        for w, b in zip(weights, bin_mids):
            Esd1.append(
                (w * stats["mnm"][b]["D"] + (1 - w) * stats["ism"][b]["D"])
                / (w * stats["mnm"][b]["pi2"] + (1 - w) * stats["ism"][b]["pi2"])
            )
            Esd2.append(
                (w * stats["mnm"][b]["D2"] + (1 - w) * stats["ism"][b]["D2"])
                / (w * stats["mnm"][b]["pi2"] + (1 - w) * stats["ism"][b]["pi2"])
            )
        return np.array(Esd1), np.array(Esd2)

    p0 = [0.1, 0.1]
    lower_bound = [0, 0]
    upper_bound = [1, 1]

    opt_params, LL = optimize_log_fmin(
        p0,
        [sd1, sd2],
        [v1, v2],
        model_func,
        verbose=0,
        lower_bound=lower_bound,
        upper_bound=upper_bound,
    )

    print("####")
    print("r:", r)
    print("opt params:", opt_params)
    print("LL:", LL)
    print("####")
    fout.write(f"{r}\t{opt_params[0]}\t{opt_params[1]}\t{LL}\n")

fout.close()
