# Simulate some number of coding regions of a given length, without the flanking
# neutral regions. Within each coding region, we have some DFE (mimicking gene
# region) with synonymous, missense, and loss-of-function mutations, which have
# mutation rates and DFEs who's defaults are from empirical studies (DFE from MSL,
# mutation rates from gnomad).


import fwdpy11
import numpy as np
import sys
import argparse
from datetime import datetime
from collections import defaultdict
import copy
import pickle
import moments.LD


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)
    sys.stderr.flush()


def current_time():
    return " [" + datetime.strftime(datetime.now(), "%Y-%m-%d %H:%M:%S") + "]"


def make_parser():
    ADHF = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser("simulate_regions.py", formatter_class=ADHF)
    parser.add_argument("--random_seed", required=True, type=int)
    optional = parser.add_argument_group("Optional")
    optional.add_argument(
        "--population_size",
        "-N",
        type=int,
        default=5000,
        help="Diploid population size, defaults to 5,000.",
    )
    optional.add_argument(
        "--recombination_rate",
        type=float,
        default=2e-8,
        help="Per-base recombination rate, defaults to 2e-8.",
    )
    optional.add_argument(
        "--mutation_rate",
        type=float,
        default=2e-8,
        help="Per-base mutation rate, defaults to 2e-8.",
    )
    optional.add_argument(
        "--sequence_length",
        default=100000,
        type=float,
        help="Length of coding region, defaults to 100,000 bp.",
    )
    optional.add_argument(
        "--burnin_factor",
        type=int,
        default=10,
        help="Burn in time, in units of 2N generations, defaults to 10.",
    )
    optional.add_argument(
        "--selection_coefficient",
        "-s",
        type=float,
        default=-0.0005,
        help="Selection coefficient for selected variants, defaults to 0.001.",
    )
    optional.add_argument(
        "--epistasis_coefficient",
        "-e",
        type=float,
        default=0,
        help="Epistasis coefficient, sAB = (sA + sB) * (1 + e), defaults to 0.",
    )
    optional.add_argument(
        "--dominance_coefficient",
        "-d",
        type=float,
        default=0.5,
        help="Dominance coefficient, defaults to 0.5 (additive).",
    )
    optional.add_argument(
        "--sample_size",
        "-n",
        type=int,
        default=100,
        help="The sample size for sampling the genotype matrix.",
    )
    optional.add_argument(
        "--sample_spacing",
        type=int,
        default=50,
        help="Spacing for how often to preserve ancient samples.",
    )
    optional.add_argument(
        "--simulation_length",
        type=int,
        default=1,
        help="Recorded simulation time, in units of 2N generations, defaults to 1.",
    )
    optional.add_argument(
        "--no_output",
        action="store_true",
        help="If flag is given, don't save anything.",
    )
    return parser


def setup_simulation(args):
    L = args.sequence_length
    # recombination
    rregions = [
        fwdpy11.PoissonInterval(0, L, args.sequence_length * args.recombination_rate)
    ]
    # selection
    U_sel = args.sequence_length * args.mutation_rate
    sregions = [
        fwdpy11.ConstantS(
            0,
            L,
            weight=U_sel,
            s=args.selection_coefficient,
            h=args.dominance_coefficient,
            label=1,
        )
    ]
    # fitness function (multiplicative vs epistasis)
    if args.epistasis_coefficient == 0:
        gvalue = fwdpy11.Multiplicative(2.0)
    pdict = {
        "gvalue": gvalue,
        "rates": (0, U_sel, None),
        "nregions": [],
        "sregions": sregions,
        "recregions": rregions,
        "demography": fwdpy11.DiscreteDemography(),
        "simlen": 2
        * args.population_size
        * (args.burnin_factor + args.simulation_length),
    }

    return pdict, L


def runsim(args):
    pdict, L = setup_simulation(args)
    pop = fwdpy11.DiploidPopulation(args.population_size, L)
    rng = fwdpy11.GSLrng(args.random_seed)

    simplification_interval = 100

    params = fwdpy11.ModelParams(**pdict)

    # sample every `sample_spacing` generations, after burnin
    ancient_sample_times = np.arange(
        2 * args.population_size * args.burnin_factor,
        2 * args.population_size * (args.burnin_factor + args.simulation_length),
        args.sample_spacing,
    )

    fwdpy11.evolvets(
        rng,
        pop,
        params,
        simplification_interval,
        recorder=fwdpy11.RandomAncientSamples(
            args.random_seed, args.sample_size, ancient_sample_times
        ),
    )
    return pop


def get_genotype_matrices(ts, args):
    sample_times = {}
    for s in ts.samples():
        if ts.node(s).time in sample_times:
            sample_times[ts.node(s).time].append(s)
        else:
            sample_times[ts.node(s).time] = [s]

    for t in sample_times.keys():
        if len(sample_times[t]) > 2 * args.sample_size:
            sample_times[t] = list(
                np.random.choice(sample_times[0], 2 * args.sample_size, replace=False)
            )
        elif len(sample_times[t]) < 2 * args.sample_size:
            raise ValueError("not enough samples!")

    Gs = {}
    positions = {}
    for t, sample_set in sample_times.items():
        ts_simp = ts.simplify(sample_set)
        positions[t] = [ts_simp.site(m.site).position for m in ts_simp.mutations()]
        Gs[t] = ts_simp.genotype_matrix()

    return positions, Gs


def get_position_distances(pos):
    # pos = list of positions, ascending
    return np.array([j - i for idx, i in enumerate(pos[:-1]) for j in pos[idx + 1 :]])


def parse_genotype_matrix(pos, G, r_bins, r):
    r_dists = r * get_position_distances(pos)
    D2, Dz, pi2, D = moments.LD.Parsing.compute_pairwise_stats(G)
    assert len(r_dists) == len(D2)
    D2_bin = np.zeros(len(r_bins) - 1)
    Dz_bin = np.zeros(len(r_bins) - 1)
    pi2_bin = np.zeros(len(r_bins) - 1)
    D_bin = np.zeros(len(r_bins) - 1)
    num_tallied = 0
    for ii, (r_l, r_r) in enumerate(zip(r_bins[:-1], r_bins[1:])):
        idx = np.where(np.logical_and(r_dists > r_l, r_dists <= r_r))[0]
        D2_bin[ii] = np.sum(D2[idx])
        Dz_bin[ii] = np.sum(Dz[idx])
        pi2_bin[ii] = np.sum(pi2[idx])
        D_bin[ii] = np.sum(D[idx])
        num_tallied += len(idx)
    assert num_tallied == len(r_dists)
    return D2_bin, Dz_bin, pi2_bin, D_bin


def parse_all_genotype_matrices(positions, Gs, r_bins, r):
    data = [
        np.zeros(len(r_bins) - 1),
        np.zeros(len(r_bins) - 1),
        np.zeros(len(r_bins) - 1),
        np.zeros(len(r_bins) - 1),
    ]
    for k in positions.keys():
        pos = positions[k]
        G = Gs[k]
        for ii, v in enumerate(parse_genotype_matrix(pos, G, r_bins, r)):
            data[ii] += v
    return data


if __name__ == "__main__":
    parser = make_parser()
    args = parser.parse_args(sys.argv[1:])

    eprint(
        current_time(),
        "Starting simulation with "
        f"{(args.burnin_factor + args.simulation_length) * 2 * args.population_size} "
        "total generations",
    )
    pop = runsim(args)
    eprint(current_time(), f"Finished simulation, at generation {pop.generation}")

    ts = pop.dump_tables_to_tskit()
    positions, Gs = get_genotype_matrices(ts, args)

    r = args.recombination_rate
    r_bins = np.array(
        [
            0,
            500 * r,
            1000 * r,
            2000 * r,
            5000 * r,
            10000 * r,
            20000 * r,
            50000 * r,
            args.sequence_length * r,
        ]
    )

    # data is all sums of [D2, Dz, pi2, D]
    data = parse_all_genotype_matrices(positions, Gs, r_bins, r)
    print(data)
