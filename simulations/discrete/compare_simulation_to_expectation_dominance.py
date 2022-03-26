import numpy as np
import argparse, sys
import moments.TwoLocus
import pickle
import matplotlib.pylab as plt


def make_parser():
    ADHF = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser("compare_simulation_to_expectation.py", formatter_class=ADHF)
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
        "--dominance_coefficient",
        "-H",
        type=float,
        default=0,
        help="The dominance coefficient at both loci.",
    )
    parser.add_argument(
        "--sample_size",
        "-n",
        type=int,
        default=50,
        help="The sample size to draw each sampling generation",
    )
    parser.add_argument(
        "--recombination_rate",
        "-r",
        type=float,
        default=0,
    )
    return parser


def get_residuals(model, data):
    """
    Return an array of Anscombe residuals, unordered.
    """
    resids = moments.Inference.Anscombe_Poisson_residual(
        data[model.mask == False], model[data.mask == False]
    )
    return -resids[resids.mask == False].data


def plot_residuals(resids, ax=None):
    """
    Plot a histogram of residuals.
    """
    if ax is None:
        fig = plt.figure()
        ax = plt.subplot(1, 1, 1)

    ax.hist(resids, bins=20, density=True)
    ax.set_xlabel("Residuals")
    ax.set_ylabel("Density")


def plot_pAB_comparison(model, data, nA=None, nB=None, ax=None, legend=False):
    if nA is None or nB is None:
        raise ValueError("must pass nA and nB")
    if ax is None:
        fig = plt.figure()
        ax = plt.subplot(1, 1, 1)

    nAB, pAB_model = moments.TwoLocus.Util.pAB(model, nA, nB)
    # variance equals expectation, assuming poisson
    pAB_error = np.sqrt(pAB_model) * 1.96
    nAB, pAB_data = moments.TwoLocus.Util.pAB(data, nA, nB)
    ax.bar(nAB - 0.2, pAB_model, width=0.35, label="E")
    ax.bar(nAB + 0.2, pAB_data, width=0.35, label="D") # , yerr=pAB_error)

    if legend:
        ax.legend(loc="upper center")
    ax.set_title(f"$n_A={nA}, n_B={nB}$")
    ax.set_ylabel("Count")
    ax.set_xlabel("$n_{AB}$")
    ax.set_xticks(nAB)


subplots_dim = (2, 4)


def plot_grid(model, data, args):
    fig = plt.figure(figsize=(8, 4))
    fig.clf()

    # flat residuals histogram
    ax1 = plt.subplot2grid(subplots_dim, (0, 0), colspan=2)
    resids = get_residuals(model, data)
    plot_residuals(resids, ax=ax1)

    nA = 1
    nB = 1
    ax2 = plt.subplot2grid(subplots_dim, (0, 2))
    plot_pAB_comparison(model, data, nA=nA, nB=nB, ax=ax2)
    ax2.set_yscale("log")

    nA = 2
    nB = 2
    ax3 = plt.subplot2grid(subplots_dim, (0, 3))
    plot_pAB_comparison(model, data, nA=nA, nB=nB, ax=ax3, legend=True)
    ax3.set_yscale("log")

    nA = 1
    nB = 2
    ax4 = plt.subplot2grid(subplots_dim, (1, 0))
    plot_pAB_comparison(model, data, nA=nA, nB=nB, ax=ax4)
    ax4.set_yscale("log")

    nA = 2
    nB = 3
    ax5 = plt.subplot2grid(subplots_dim, (1, 1))
    plot_pAB_comparison(model, data, nA=nA, nB=nB, ax=ax5)
    ax5.set_yscale("log")

    nA = 3
    nB = 3
    ax6 = plt.subplot2grid(subplots_dim, (1, 2))
    plot_pAB_comparison(model, data, nA=nA, nB=nB, ax=ax6)
    ax6.set_yscale("log")

    nA = 4
    nB = 5
    ax7 = plt.subplot2grid(subplots_dim, (1, 3))
    plot_pAB_comparison(model, data, nA=nA, nB=nB, ax=ax7)
    ax7.set_yscale("log")

    fig.suptitle(rf"$N={args.population_size}, n={args.sample_size}, s=-{args.selection_coefficient}, h={args.dominance_coefficient}, r={args.recombination_rate}$")
    fig.tight_layout()


if __name__ == "__main__":
    parser = make_parser()
    args = parser.parse_args(sys.argv[1:])

    Ne = args.population_size
    n = args.sample_size
    r = args.recombination_rate
    s = args.selection_coefficient
    h = args.dominance_coefficient

    try:
        fname = f"outputs/Ne_{Ne}_n_{n}_r_{r}_s_{s}_h_{h}.bp"
        all_data = pickle.load(
            open(fname, "rb")
        )
    except IOError:
        try:
            if e == 0.0:
                e = 0
            fname = f"outputs/Ne_{Ne}_n_{n}_r_{r}_s_{s}_h_{h}.bp"
            all_data = pickle.load(
                open(fname, "rb")
            )
        except IOError:
            raise IOError(f"data does not exist for file name: {fname}")

    data = all_data["simulation"]

    n_proj = 80
    fname = f"expectations/Ne_{Ne}_from_{n_proj}_n_{n}_r_{r}_s_{s}_h_{h}.bp"
    try:
        model = moments.TwoLocus.TLSpectrum.from_file(fname)
    except IOError:
        theta = 4 * Ne * all_data["args"].mutation_rate
        rho = 4 * Ne * r
        gamma = - 2 * Ne * s
        sel_params_general = moments.TwoLocus.Util.simple_dominance(gamma, h=h)

        model = moments.TwoLocus.Demographics.equilibrium(
            n_proj,
            rho=rho,
            theta=theta,
            sel_params_general=sel_params_general
        ) * all_data["args"].num_replicates

        model = model.project(n)
        model.to_file(fname)

    data.mask_fixed()
    model.mask_fixed()

    plot_grid(model, data, all_data["args"])
    plt.savefig(f"plots/comp_Ne_{Ne}_n_{n}_r_{r}_s_{s}_h_{h}.pdf")
    #plt.show()
