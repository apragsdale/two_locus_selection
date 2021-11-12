# For a given population size, selection coefficient, dominance coefficient or
# epistasis coefficient, and sample size, run some number of


import moments.TwoLocus
import numpy as np
import sys, os
import pickle


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
        default=5000,
        help="Diploid population size, defaults to 5,000.",
    )
    parser.add_argument(
        "--selection_coefficient",
        "-s",
        type=float,
        default=0,
        help="gamma = -2 * Ne * s. Note that we set the abs value",
    )
    parser.add_argument(
        "--epistasis_coefficient",
        "-e",
        type=float,
        default=0,
        help="The epistasis coefficient, so that s_{AB} = s_A + s_B + ?",
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
        default=1e-3,
    )
    parser.add_argument(
        "--recombination_rate",
        "-r",
        type=float,
        default=0,
    )
    parser.add_argument(
        "--num_replicates",
        type=int,
        default=100000,
    )
    parser.add_argument(
        "--spacing",
        type=int,
        default=100,
    )
    parser.add_argument(
        "--out",
        "-o",
        action="store_true",
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


def sel_rec(X, sA, sB, sAB, r):
    x = X / 2 / Ne
    sel_fac = 1 / (1 + sAB * x[:, 0] + sA * x[:, 1] + sB * x[:, 2])
    xAB = x[:, 0] * (1 + sAB) * sel_fac
    xAb = x[:, 1] * (1 + sA) * sel_fac
    xaB = x[:, 2] * (1 + sB) * sel_fac
    xab = x[:, 3] * sel_fac
    # recombination
    pAB = xAB ** 2 + xAB * xAb + xAB * xaB + (1 - r) * xAB * xab + r * xAb * xaB
    pAb = xAB * xAb + r * xAB * xab + xAb ** 2 + (1 - r) * xAb * xaB + xAb * xab
    paB = xAB * xaB + r * xAB * xab + (1 - r) * xaB * xAb + xaB ** 2 + xaB * xab
    pab = (1 - r) * xAB * xab + r * xAb * xaB + xAb * xab + xaB * xab + xab ** 2
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


def run_sim(Ne, n, theta, sAB, sA, sB, args):
    F = np.zeros((n + 1, n + 1, n + 1))
    FA = np.zeros(n + 1)
    FB = np.zeros(n + 1)

    X = np.empty((0, 4), dtype=np.int)
    XA = np.empty((0,), dtype=np.int)
    XB = np.empty((0,), dtype=np.int)

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

        ## single locus evolution
        pA = XA * (1 + sA) / (2 * Ne + sA * XA)
        XA = np.random.binomial(2 * Ne, pA)
        pB = XB * (1 + sB) / (2 * Ne + sB * XB)
        XB = np.random.binomial(2 * Ne, pB)

        ## two locus evolution
        # selection and recombination
        p = sel_rec(X, sA, sB, sAB, r)

        # draw new generation
        X = draw_generation(X, p, Ne)

        # remove fixed pairs
        X = remove_fixed(X)
        XA = XA.compress(np.logical_and(XA != 0, XA != 2 * Ne))
        XB = XB.compress(np.logical_and(XB != 0, XB != 2 * Ne))

        if (gen + 1) % Ne == 0:
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
    eprint("Input arguments:", args)

    n = args.sample_size
    Ne = args.population_size
    rho = 4 * Ne * args.recombination_rate
    theta = 4 * Ne * args.mutation_rate

    sA = sB = -args.selection_coefficient
    sAB = (sA + sB) * (1 + args.epistasis_coefficient)

    eprint("Scaled parameters:")
    eprint("gamma:", -2 * Ne * args.selection_coefficient)
    eprint("rho:", 4 * Ne * args.recombination_rate)
    eprint("theta:", 4 * Ne * args.mutation_rate)
    F = run_sim(Ne, n, theta, sAB, sA, sB, args)

    #sel_params = moments.TwoLocus.Util.additive_epistasis(
    #    sA, epsilon=args.epistasis_coefficient, Ne=Ne
    #)

    #E_F = moments.TwoLocus.Demographics.equilibrium(
    #    n, rho, sel_params=sel_params, theta=theta
    #) * args.num_replicates

    if args.out is True:
        with open(
            f"outputs/Ne_{Ne}_n_{args.sample_size}_s_{args.selection_coefficient}_e_{args.epistasis_coefficient}.bp",
            "wb+",
        ) as fout:
            # pickle.dump({"args": args, "simulation": F, "expectation": E_F}, fout)
            pickle.dump({"args": args, "simulation": F}, fout)
