import moments
import pickle
import numpy as np
import matplotlib.pylab as plt

populations = {
    "Africa": ["YRI"],  # , "GWD", "LWK", "MSL", "YRI"],
    "Europe": ["CEU"],  # , "FIN", "GBR", "IBS", "TSI"],
    "East Asia": ["CHB"],  # , "CHB", "CHS", "JPT", "KHV"],
}

import bokeh.palettes

colors = bokeh.palettes.Paired[8]


def plot_mis_syn(ax, mis, syn, title=None, xlabel=False, ylabel=False, legend=False):
    mis = mis.fold()
    syn = syn.fold()

    ax.plot(mis / mis.S(), ".-", color=colors[6], ms=2, lw=1, label="Missense (all)")
    ax.plot(syn / syn.S(), ".--", color=colors[0], ms=2, lw=1, label="Synonymous (all)")

    mis.mask[: int(mis.sample_sizes * 0.1)] = True
    syn.mask[: int(mis.sample_sizes * 0.1)] = True
    ax.plot(mis / mis.S(), ".-", color=colors[7], ms=2, lw=1, label="Missense (>10%)")
    ax.plot(
        syn / syn.S(), ".--", color=colors[1], ms=2, lw=1, label="Synonymous (>10%)"
    )

    if legend:
        ax.legend(fontsize=6)
    ax.set_yscale("log")
    if ylabel:
        ax.set_ylabel("Proportion")
    if xlabel:
        ax.set_xlabel("Minor allele frequency")
    ax.set_title(title)


fig = plt.figure(figsize=(8, 8))

i = 0

for c, ps in populations.items():
    for p in ps:
        fname = f"../analysis/parsed_data/{p}.frequency_spectra.bp"
        data = pickle.load(open(fname, "rb"))

        i += 1
        ax = plt.subplot(3, 3, i)
        if i > 6:
            xlabel = True
        else:
            xlabel = False
        if i == 1:
            legend = True
        else:
            legend = False

        # all
        mis = moments.Spectrum(data["missense"]["all"])
        syn = moments.Spectrum(data["synonymous"]["all"])
        plot_mis_syn(
            ax,
            mis,
            syn,
            title=f"{p}, all mutations",
            ylabel=True,
            xlabel=xlabel,
            legend=legend,
        )

        i += 1
        # within
        ax = plt.subplot(3, 3, i)
        mis = moments.Spectrum(data["missense"]["in"])
        syn = moments.Spectrum(data["synonymous"]["in"])
        plot_mis_syn(
            ax,
            mis,
            syn,
            title=f"{p}, within domains",
            xlabel=xlabel,
        )

        i += 1
        # outside
        ax = plt.subplot(3, 3, i)
        mis = moments.Spectrum(data["missense"]["out"])
        syn = moments.Spectrum(data["synonymous"]["out"])
        plot_mis_syn(
            ax,
            mis,
            syn,
            title=f"{p}, outside domains",
            xlabel=xlabel,
        )

fig.tight_layout()
plt.savefig("SFS_proportions.pdf", dpi=300)
plt.show()
