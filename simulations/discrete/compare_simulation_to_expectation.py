import numpy as np
import moments.TwoLocus
import pickle
import matplotlib.pylab as plt


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


def plot_pAB_comparison(model, data, nA=None, nB=None, ax=None):
    if nA is None or nB is None:
        raise ValueError("must pass nA and nB")
    if ax is None:
        fig = plt.figure()
        ax = plt.subplot(1, 1, 1)

    nAB, pAB_model = moments.TwoLocus.Util.pAB(model, nA, nB)
    nAB, pAB_data = moments.TwoLocus.Util.pAB(data, nA, nB)
    ax.bar(nAB - 0.2, pAB_model / pAB_model.sum(), width=0.35, label="E[pAB]")
    ax.bar(nAB + 0.2, pAB_data / pAB_data.sum(), width=0.35, label="O[pAB]")

    ax.set_title(f"$n_A={nA}, n_B={nB}$")
    ax.set_ylabel("Probability")
    ax.set_xlabel("$n_{AB}$")
    ax.set_xticks(nAB)


subplots_dim = (2, 4)


def plot_grid(model, data, args):
    fig = plt.figure()

    # flat residuals histogram
    ax1 = plt.subplot2grid(subplots_dim, (0, 0), colspan=2)
    resids = get_residuals(model, data)
    plot_residuals(resids, ax=ax1)

    # nA = 1, nB = 1
    ax2 = plt.subplot2grid(subplots_dim, (0, 2))
    plot_pAB_comparison(model, data, nA=1, nB=1, ax=ax2)

    # nA = 2, nB = 3
    ax3 = plt.subplot2grid(subplots_dim, (0, 3))
    plot_pAB_comparison(model, data, nA=2, nB=3, ax=ax3)

    # nA = 10, nB = 15
    ax4 = plt.subplot2grid(subplots_dim, (1, 0), colspan=2)
    plot_pAB_comparison(model, data, nA=4, nB=4, ax=ax4)

    # nA = 8, nB = 30
    ax4 = plt.subplot2grid(subplots_dim, (1, 2), colspan=2)
    plot_pAB_comparison(model, data, nA=8, nB=8, ax=ax4)

    fig.suptitle("[[details of sim]]")
    fig.tight_layout()


if __name__ == "__main__":
    all_data = pickle.load(
        open("outputs/Ne_100_n_10_r_0.0025_s_0_e_0.bp", "rb")
    )
    data = all_data["simulation"]
    args = all_data["args"]

    Ne = args.population_size
    theta = 4 * Ne * args.mutation_rate
    rho = 4 * Ne * args.recombination_rate
    gamma = - 2 * Ne * args.selection_coefficient
    sel_params = [2 * gamma * (1 + args.epistasis_coefficient), gamma, gamma]

    model = moments.TwoLocus.Demographics.equilibrium(
        args.sample_size,
        rho=rho,
        theta=theta,
        sel_params=sel_params
    ) * args.num_replicates

    data.mask_fixed()
    model.mask_fixed()

    plot_grid(model, data, all_data["args"])
    plt.show()
