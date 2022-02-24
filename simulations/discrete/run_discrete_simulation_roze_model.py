import moments.TwoLocus
import numpy as np
import sys, os
import pickle
import gzip
import time

import argparse
from datetime import datetime


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)
    sys.stderr.flush()


def current_time():
    return " [" + datetime.strftime(datetime.now(), "%Y-%m-%d %H:%M:%S") + "]"


def make_parser():
    ADHF = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser("run_discrete_comparison.py", formatter_class=ADHF)
    parser.add_argument(
        "--population_size",
        "-N",
        type=int,
        default=1000,
        help="Diploid population size, defaults to 5,000.",
    )
    parser.add_argument(
        "--compound_coefficient",
        "-s",
        type=float,
        default=0.01,
    )
    parser.add_argument(
        "--dominance_coefficient",
        "-d",
        type=float,
        default=0.1,
        help="The dominance coefficient, so that s_{AA} = 2 * d * s_A",
    )
    parser.add_argument(
        "--sample_size",
        "-n",
        type=int,
        default=50,
        help="The sample size to draw each sampling generation",
    )
    parser.add_argument(
        "--mutation_rate",
        "-u",
        type=float,
        default=1e-4,
    )
    parser.add_argument(
        "--recombination_rate",
        "-r",
        type=float,
        default=0.0001,
    )
    parser.add_argument(
        "--num_replicates",
        type=int,
        default=500000,
    )
    parser.add_argument(
        "--spacing",
        type=int,
        default=50,
    )
    parser.add_argument(
        "--out",
        "-o",
        action="store_true",
    )
    parser.add_argument(
        "--index",
        "-i",
        default=None,
        type=int,
    )
    return parser


def draw_mutation_pairs_discrete(U, XA, XB):
    if len(XA) > 0:
        nA = np.random.choice(XA, size=np.random.poisson(len(XA) * U / 2.0))
        nAB = (np.random.rand(len(nA)) < nA / 2 / Ne).astype("int")
        nAb = nA - nAB
        naB = np.ones(len(nA), dtype=np.int) - nAB
        nab = 2 * Ne - nAB - nAb - naB
        new_pairs_A = np.array([nAB, nAb, naB, nab]).T
    else:
        new_pairs_A = np.empty((0, 4), dtype=np.int)
    if len(XB) > 0:
        nB = np.random.choice(XB, size=np.random.poisson(len(XB) * U / 2.0))
        nAB = (np.random.rand(len(nB)) < nB / 2 / Ne).astype("int")
        naB = nB - nAB
        nAb = np.ones(len(nB), dtype=np.int) - nAB
        nab = 2 * Ne - nAB - nAb - naB
        new_pairs_B = np.array([nAB, nAb, naB, nab]).T
    else:
        new_pairs_B = np.empty((0, 4), dtype=np.int)
    new_pairs = np.concatenate((new_pairs_A, new_pairs_B))
    return new_pairs


def sel_rec(X, S, r):
    # S from moments.TwoLocus.Util.simple_dominance(s, h), so
    # S = [sABAB, sABAb, sABaB, sABab, sAbAb, sAbaB, sAbab, saBaB, saBab]
    # note that this could be easily extended to account for gene-based dominance
    [sABAB, sABAb, sABaB, sABab, sAbAb, sAbaB, sAbab, saBaB, saBab] = S
    x = X / 2 / Ne
    # make diploids
    xABAB = x[:, 0] ** 2
    xABAb = 2 * x[:, 0] * x[:, 1]
    xABaB = 2 * x[:, 0] * x[:, 2]
    xABab = 2 * x[:, 0] * x[:, 3]
    xAbAb = x[:, 1] ** 2
    xAbaB = 2 * x[:, 1] * x[:, 2]
    xAbab = 2 * x[:, 1] * x[:, 3]
    xaBaB = x[:, 2] ** 2
    xaBab = 2 * x[:, 2] * x[:, 3]
    xabab = x[:, 3] ** 2
    assert np.allclose(
        xABAB + xABAb + xABaB + xABab + xAbAb + xAbaB + xAbab + xaBaB + xaBab + xabab,
        1
    )
    # select on diploids
    sel_fac = 1 / (
        1 + sABAB * xABAB + sABAb * xABAb + sABaB * xABaB + sABab * xABab
        + sAbAb * xAbAb + sAbaB * xAbaB + sAbab * xAbab
        + saBaB * xaBaB + saBab * xaBab
    )
    yABAB = xABAB * (1 + sABAB) * sel_fac
    yABAb = xABAb * (1 + sABAb) * sel_fac
    yABaB = xABaB * (1 + sABaB) * sel_fac
    yABab = xABab * (1 + sABab) * sel_fac
    yAbAb = xAbAb * (1 + sAbAb) * sel_fac
    yAbaB = xAbaB * (1 + sAbaB) * sel_fac
    yAbab = xAbab * (1 + sAbab) * sel_fac
    yaBaB = xaBaB * (1 + saBaB) * sel_fac
    yaBab = xaBab * (1 + saBab) * sel_fac
    yabab = xabab * sel_fac
    assert np.allclose(
        yABAB + yABAb + yABaB + yABab + yAbAb + yAbaB + yAbab + yaBaB + yaBab + yabab,
        1
    )
    # recombination
    pAB = yABAB + yABAb / 2 + yABaB / 2 + (1 - r) * yABab / 2 + r * yAbaB / 2
    pAb = yABAb / 2 + r * yABab / 2 + yAbAb + (1 - r) * yAbaB / 2 + yAbab / 2
    paB = yABaB / 2 + r * yABab / 2 + (1 - r) * yAbaB / 2 + yaBaB + yaBab / 2
    pab = (1 - r) * yABab / 2 + r * yAbaB / 2 + yAbab / 2 + yaBab / 2 + yabab
    # check sums to one
    sums = pAB + pAb + paB + pab
    assert np.allclose(sums, 1)
    return np.array([pAB, pAb, paB, pab]).T


def draw_generation(X, p, Ne):
    for ii, p_ii in enumerate(p):
        X[ii] = np.random.multinomial(2 * Ne, p_ii)
    return X


def remove_fixed(X):
    nA = X[:, 0] + X[:, 1]
    nB = X[:, 0] + X[:, 2]
    not_fixed = 1 - np.logical_or(
        np.logical_or(nA == 0, nA == 2 * Ne), np.logical_or(nB == 0, nB == 2 * Ne)
    )
    X = X.compress(not_fixed, axis=0)
    return X


#generator = np.random.Generator(np.random.PCG64(np.random.randint(1000000)))


def sample(X, ns, F):
    for ii, x in enumerate(X):
        [nAB, nAb, naB, nab] = np.random.multinomial(ns, x / 2 / Ne)
        # [nAB, nAb, naB, nab] = generator.multivariate_hypergeometric(x, ns)
        F[nAB, nAb, naB] += 1
    return F


def sample_single(X, Ne, ns, F):
    # i = np.random.binomial(ns, X / 2 / Ne)
    i = np.random.hypergeometric(X, 2 * Ne - X, ns)
    inds, counts = np.unique(i, return_counts=True)
    F[inds] += counts
    return F


def single_locus_sel(x, s, h):
    xA = x ** 2 * (1 + s) ** 2 + 2 * x * (1 - x) * (1 + 2 * h * s) / 2
    xa = 2 * x * (1 - x) * (1 + 2 * h * s) / 2 + (1 - x) ** 2
    return xA / (xA + xa)

def run_sim(Ne, n, theta, s, h, args):
    F = np.zeros((n + 1, n + 1, n + 1))
    FA = np.zeros(n + 1)
    FB = np.zeros(n + 1)

    X = np.empty((0, 4), dtype=np.int)
    XA = np.empty((0,), dtype=np.int)
    XB = np.empty((0,), dtype=np.int)

    S = moments.TwoLocus.Util.simple_dominance(s, h=h)
    r = args.recombination_rate

    burnin_gens = 160 * Ne
    total_gens = burnin_gens + args.num_replicates * args.spacing
    eprint(current_time(), "starting simulation")
    for gen in range(total_gens):
        # draw single site background mutations
        XA = np.concatenate((XA, np.ones(np.random.poisson(theta / 2.0), dtype=np.int)))
        XB = np.concatenate((XB, np.ones(np.random.poisson(theta / 2.0), dtype=np.int)))

        # draw new pairs of mutations
        new_pairs = draw_mutation_pairs_discrete(theta, XA, XB)
        X = np.concatenate((X, new_pairs))

        # single locus evolution
        pA = single_locus_sel(XA / 2 / Ne, s, h)
        XA = np.random.binomial(2 * Ne, pA)
        pB = single_locus_sel(XB / 2 / Ne, s, h)
        XB = np.random.binomial(2 * Ne, pB)

        ## two locus evolution
        # selection and recombination
        p = sel_rec(X, S, r)

        # draw new generation
        X = draw_generation(X, p, Ne)

        # remove fixed pairs
        X = remove_fixed(X)
        XA = XA.compress(np.logical_and(XA != 0, XA != 2 * Ne))
        XB = XB.compress(np.logical_and(XB != 0, XB != 2 * Ne))

        if (gen + 1) % (10 * Ne) == 0:
            eprint(current_time(), "at generation", gen + 1, "of", total_gens)
        if gen + 1 > burnin_gens and (gen + 1) % args.spacing == 0:
            F = sample(X, n, F)
            FA = sample_single(XA, Ne, n, FA)
            FB = sample_single(XB, Ne, n, FB)

    F = moments.TwoLocus.TLSpectrum(F)
    F[0, :, 0] = FA
    F[0, 0, :] = FB
    F.mask[0, 0, 0] = F.mask[-1, 0, 0] = F.mask[0, -1, 0] = F.mask[0, 0, -1] = True
    F.mask_infeasible()
    return F


if __name__ == "__main__":
    parser = make_parser()
    args = parser.parse_args(sys.argv[1:])
    
    if args.index is not None:
        i = args.index
        # get sh, dominance, and r from external file
        arg_options = []
        with open("inputs.txt", "r") as fin:
            for line in fin:
                arg_options.append([float(x) for x in line.split()])
        args.compound_coefficient = arg_options[i - 1][0]
        args.dominance_coefficient = arg_options[i - 1][1]
        args.recombination_rate = arg_options[i - 1][2]
    
    eprint("Input arguments:", args)

    n = args.sample_size
    Ne = args.population_size
    rho = 4 * Ne * args.recombination_rate
    theta = 4 * Ne * args.mutation_rate

    sh = args.compound_coefficient
    h = args.dominance_coefficient
    s = -sh / h

    eprint("Scaled parameters:")
    eprint("gamma:", 2 * Ne * s)
    eprint("dominance:", h)
    eprint("rho:", rho)
    eprint("theta:", theta)
    
    F = run_sim(Ne, n, theta, s, h, args)

    eprint("Finished simulation")

    if args.out is True:
        t = int(time.time())
        fname = f"outputs/Ne_{Ne}_n_{args.sample_size}_r_{args.recombination_rate}_sh_{args.compound_coefficient}_h_{args.dominance_coefficient:0.2f}_time_{t}.bp.gz"
        with gzip.open(fname, "wb+") as fout:
            # pickle.dump({"args": args, "simulation": F, "expectation": E_F}, fout)
            pickle.dump({"args": args, "simulation": F}, fout)
