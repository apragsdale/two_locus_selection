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
import matplotlib.pylab as plt
import msprime, tskit
import gzip


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
        default=1000,
        help="Diploid population size, defaults to 5,000.",
    )
    optional.add_argument(
        "--recombination_rate",
        type=float,
        default=1e-8,
        help="Per-base recombination rate, defaults to 2e-8.",
    )
    optional.add_argument(
        "--mutation_rate",
        type=float,
        default=1e-8,
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
        default=25,
        help="The sample size for sampling the genotype matrix.",
    )
    optional.add_argument(
        "--sample_spacing",
        type=int,
        default=100,
        help="Spacing for how often to preserve ancient samples.",
    )
    optional.add_argument(
        "--simulation_length",
        type=int,
        default=100,
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
            h=2 * args.dominance_coefficient,
            label=1,
        )
    ]
    # fitness function (multiplicative vs epistasis)
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
    ## TODO: also get genotype matrices for neutral mutations at same mutation rate
    ## dropping them in after clearing the selected mutations
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
    init_num_muts = ts.num_mutations

    Gs = {}
    Gs_neu = {}
    positions = {}
    positions_neu = {}
    for t, sample_set in sample_times.items():
        ts_simp = ts.simplify(sample_set)
        positions[t] = [ts_simp.site(m.site).position for m in ts_simp.mutations()]
        Gs[t] = ts_simp.genotype_matrix()

        # clear mutations
        ts_simp = ts_simp.delete_sites(range(ts_simp.num_sites))
        assert ts_simp.num_mutations == 0
        # add neutral mutations
        ts_simp = msprime.sim_mutations(
            ts_simp, rate=args.mutation_rate, discrete_genome=False
        )
        # get Gs_neutral
        positions_neu[t] = [ts_simp.site(m.site).position for m in ts_simp.mutations()]
        Gs_neu[t] = ts_simp.genotype_matrix()
    assert ts.num_mutations == init_num_muts
    return positions, Gs, positions_neu, Gs_neu


def get_position_distances(pos):
    # pos = list of positions, ascending
    return np.array([j - i for idx, i in enumerate(pos[:-1]) for j in pos[idx + 1 :]])


def parse_genotype_matrix(pos, G, r_bins, r):
    # r_dists = r * get_position_distances(pos)
    # D2, Dz, pi2, D = moments.LD.Parsing.compute_pairwise_stats(G)
    # assert len(r_dists) == len(D2)
    # D2_bin = np.zeros(len(r_bins) - 1)
    # Dz_bin = np.zeros(len(r_bins) - 1)
    # pi2_bin = np.zeros(len(r_bins) - 1)
    # D_bin = np.zeros(len(r_bins) - 1)
    # num_tallied = 0
    # for ii, (r_l, r_r) in enumerate(zip(r_bins[:-1], r_bins[1:])):
    # idx = np.where(np.logical_and(r_dists > r_l, r_dists <= r_r))[0]
    # D2_bin[ii] = np.sum(D2[idx])
    # Dz_bin[ii] = np.sum(Dz[idx])
    # pi2_bin[ii] = np.sum(pi2[idx])
    # D_bin[ii] = np.sum(D[idx])
    # num_tallied += len(idx)
    # assert num_tallied == len(r_dists)
    psi = [
        np.zeros((G.shape[1] + 1, G.shape[1] + 1, G.shape[1] + 1))
        for _ in range(len(r_bins) - 1)
    ]
    for ii, (pos_l, gs_l) in enumerate(zip(pos, G)):
        for jj, (pos_r, gs_r) in enumerate(zip(pos, G)):
            if jj <= ii:
                continue
            r_dist = r * (pos_r - pos_l)
            bin_idx = np.where(
                np.logical_and(r_dist < r_bins[1:], r_dist >= r_bins[:-1])
            )[0][0]
            gs = list(zip(gs_l, gs_r))
            nAB, nAb, naB = gs.count((1, 1)), gs.count((1, 0)), gs.count((0, 1))
            psi[bin_idx][nAB, nAb, naB] += 1

    assert np.sum(psi) == len(pos) * (len(pos) - 1) / 2

    fs = np.zeros(G.shape[1] + 1)
    freqs, counts = np.unique(G.sum(axis=1), return_counts=True)
    fs[freqs] += counts
    assert np.sum(fs) == len(pos)
    # return D2_bin, Dz_bin, pi2_bin, D_bin, fs
    return psi + [fs]


def parse_all_genotype_matrices(positions, Gs, r_bins, r):
    n = Gs[0].shape[1]
    data = [
        np.zeros((n + 1, n + 1, n + 1)) for _ in range(len(r_bins) - 1)  # 2-loc fs
    ] + [np.zeros(n + 1)]

    for k in positions.keys():
        pos = positions[k]
        G = Gs[k]
        for ii, v in enumerate(parse_genotype_matrix(pos, G, r_bins, r)):
            data[ii] += v
    return data


def plot_spectra(fs_sel, fs_neu, args):
    fs_sel = moments.Spectrum(fs_sel)
    fs_neu = moments.Spectrum(fs_neu)

    # expected sfs
    reps = int(
        args.simulation_length * args.population_size * 2 / args.sample_spacing + 1
    )
    theta = reps * 4 * args.population_size * args.mutation_rate * args.sequence_length
    gamma = 2 * args.population_size * args.selection_coefficient
    h = args.dominance_coefficient
    model_sel = (
        moments.LinearSystem_1D.steady_state_1D(2 * args.sample_size, gamma=gamma, h=h)
        * theta
    )
    model_sel = moments.Spectrum(model_sel)
    model_neu = (
        moments.Spectrum(moments.LinearSystem_1D.steady_state_1D(2 * args.sample_size))
        * theta
    )

    fig = plt.figure(figsize=(6, 4))
    fig.clf()
    ax1 = plt.subplot(1, 2, 1)
    ax1.plot(fs_sel, "o", color="red", ms=3, lw=0.5, fillstyle="none", label="Data")
    ax1.plot(model_sel, ".--", color="k", ms=2, lw=0.5, label="Model")
    ax1.legend(fontsize=6)
    ax1.set_ylabel("Count")
    ax1.set_title("Selected mutations")

    top = np.max([fs_sel.max(), model_sel.max(), fs_neu.max(), model_neu.max()])
    top *= 1.2
    bottom = np.max(
        [1, np.min([fs_sel.min(), model_sel.min(), fs_neu.min(), model_neu.min()])]
    )
    bottom *= 0.1
    ax1.set_ylim(bottom=bottom, top=top)

    ax2 = plt.subplot(1, 2, 2, sharey=ax1)
    ax2.plot(fs_neu, "o", color="red", ms=3, lw=0.5, fillstyle="none", label="Data")
    ax2.plot(model_neu, ".--", color="k", ms=2, lw=0.5, label="Model")
    ax2.set_ylabel("Count")
    ax2.set_xlabel("Derived allele count")
    ax2.set_title("Neutral mutations")
    ax1.set_yscale("log")
    ax2.set_yscale("log")
    fig.tight_layout()
    return fig


if __name__ == "__main__":
    parser = make_parser()
    args = parser.parse_args(sys.argv[1:])

    eprint(current_time(), "Random seed:", args.random_seed)
    eprint(
        current_time(),
        "Parameters are "
        f"gamma = {2 * args.population_size * args.selection_coefficient}, "
        f"h = {args.dominance_coefficient}, with N = {args.population_size}",
    )
    eprint(
        current_time(),
        "Starting simulation with "
        f"{(args.burnin_factor + args.simulation_length) * 2 * args.population_size} "
        "total generations",
    )
    pop = runsim(args)
    eprint(current_time(), f"Finished simulation, at generation {pop.generation}")

    ts = pop.dump_tables_to_tskit()
    positions_sel, Gs_sel, positions_neu, Gs_neu = get_genotype_matrices(ts, args)

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

    eprint(current_time(), f"Finished getting genotype matrices, now parsing")
    # data is all sums of [D2, Dz, pi2, D]
    data_sel = parse_all_genotype_matrices(positions_sel, Gs_sel, r_bins, r)
    data_neu = parse_all_genotype_matrices(positions_neu, Gs_neu, r_bins, r)
    eprint(current_time(), "Collected", len(Gs_sel), "time points")

    gamma = 2 * args.population_size * args.selection_coefficient
    h = args.dominance_coefficient
    fname = f"data-gamma_{gamma}-h_{h}-seed_{args.random_seed}"
    pickle.dump(
        {"args": args, "data_sel": data_sel, "data_neu": data_neu},
        gzip.open("outputs/" + fname + ".bp.gz", "wb+"),
    )

    # plot frequency spectra
    #fig = plot_spectra(data_sel[-1], data_neu[-1], args)
    #plt.savefig("outputs/" + fname + ".pdf", dpi=300)
    # fig.show()
