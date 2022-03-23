import numpy as np
import itertools
import gzip
import moments.LD
import sys
import pickle
from collections import defaultdict


lof_annotations = [
    "frameshift_variant",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "start_lost",
    "stop_gained",
    "stop_lost",
    "transcript_ablation",
]

import tarfile

try:
    print("already loaded ancestral seqs", len(ancestral_seqs))
except NameError:
    ancestral_seqs = {}
    print("reading ancestral sequences")

    tar = tarfile.open(
        "/home/aaronragsdale/Data/1000G/human_ancestor/human_ancestor_GRCh37_e59.tar.bz2",
        mode="r:bz2",
    )
    for member in tar.getmembers():
        if member.isfile() is False:
            # print(member)
            continue
        if member.get_info()["name"].endswith("fa"):
            f = tar.extractfile(member)
            ## can't get SeqIO.parse to read in the extracted file
            lines = f.readlines()
            chrom = lines[0].decode().split(":")[2]
            print("found chrom", chrom)
            ancestral_seqs[chrom] = "".join([l.decode().strip() for l in lines[1:]])
            # for record in SeqIO.parse(x, "fasta"):
            #    chrom = record.id.split(":")[2]
            #    print(chrom)
            #    ancestral_seqs[chrom] = record.seq
    tar.close()


def get_annotations(annot_file):
    annots = {}  # pos: alt: (gene, csq)
    for line in open(annot_file, "r"):
        chrom, pos, ref, alt, gene, csq = line.split()
        annots.setdefault(pos, {})
        annots[pos][alt] = (gene, csq)
    return annots


def flip_gts(gts):
    gts_out = []
    for gt in gts:
        if gt == "0|0":
            gts_out.append("1|1")
        elif gt == "0|1":
            gts_out.append("1|0")
        elif gt == "1|0":
            gts_out.append("0|1")
        elif gt == "1|1":
            gts_out.append("0|0")
        else:
            raise ValueError("unexpected genotype", gt)
    return gts_out


def to_genotypes(gts):
    g = np.zeros(len(gts))
    for ii, gg in enumerate(gts):
        g[ii] = int(gg[0]) + int(gg[2])
    return g


SNPs = ["A", "C", "G", "T"]


def load_genotype_matrix(chrom, pop, haplotypes=True):
    """
    Returns positions (Lx1), annotations (Lx1), genes (Lx1), and G (Lxn)
    where L is number of loci, n number of samples
    """
    annot_file = f"./data/supp/variants/coding_variant_consequences.chr{chrom}.txt"
    annot_dict = get_annotations(annot_file)
    vcf = f"./data/vcfs/coding.{pop}.chr{chrom}.vcf.gz"
    positions = []
    annotations = []
    genes = []
    G = []
    with gzip.open(vcf, "rb") as fin:
        for line in fin:
            l = line.decode()
            if l.startswith("#"):
                continue
            pos = l.split()[1]
            ref = l.split()[3]
            alt = l.split()[4]
            if ref not in SNPs or alt not in SNPs:
                continue

            try:
                AA = l.split()[7].split(";AA=")[1].split("|")[0]
            except IndexError:
                # AA not found in info
                AA = "-"
            if AA not in SNPs:
                continue
            # get from ancestral fastas
            AA_fasta = ancestral_seqs[str(chrom)][int(pos) - 1]
            # only keep high confidence ancestral alleles
            # AA_fasta = AA_fasta.upper()
            if AA_fasta != AA:
                continue
            gts = l.split()[9:]
            if AA != ref:
                if AA == alt:
                    # flip SNP
                    gts = flip_gts(gts)
                else:
                    # doesn't match ref or alt
                    # print(f"AA = {AA}, ref = {ref}, alt = {alt}, pos = {pos}")
                    continue
            if np.all([g == "1|1" for g in gts]) or np.all([g == "0|0" for g in gts]):
                # fixed for ref or alt
                continue
            gene, csq = annot_dict[pos][alt]
            positions.append(int(pos))
            annotations.append(csq)
            genes.append(gene)
            G.append(to_genotypes(gts))

    return np.array(positions), np.array(annotations), np.array(genes), np.array(G)


def subset_annotation(positions, annotations, genes, G, annot="synonymous"):
    """
    annot options are synonymous, missense, and loss_of_function
    """
    if annot not in ["synonymous", "missense", "loss_of_function"]:
        raise ValueError("pick one of synonymous, missense, and loss_of_function")
    to_keep = np.zeros(len(positions))
    for ii, pos_annot in enumerate(annotations):
        is_lof = False
        is_mis = False
        is_syn = False
        for lof_annot in lof_annotations:
            if pos_annot.startswith(lof_annot):
                is_lof = True
        if not is_lof:
            if pos_annot.startswith("missense"):
                is_mis = True
        if not is_lof and not is_mis:
            if pos_annot.startswith("synonymous"):
                is_syn = True
        if annot == "loss_of_function" and is_lof:
            to_keep[ii] = 1
        elif annot == "missense" and is_mis:
            to_keep[ii] = 1
        elif annot == "synonymous" and is_syn:
            to_keep[ii] = 1

    new_positions = positions.compress(to_keep)
    new_genes = genes.compress(to_keep)
    new_G = G.compress(to_keep, axis=0)
    return new_positions, new_genes, new_G


def get_two_locus_genotype_counts_within_domains(genes, G, pos, dom):
    """
    For each gene, return genotype counts for each configuration seen.
    """
    genes_set = set(genes)
    genes_genotype_counts = defaultdict(list)
    genes_positions = defaultdict(list)
    for gene in genes_set:
        if gene not in dom.keys():
            continue

        to_keep = genes == gene
        G_gene = G.compress(to_keep, axis=0)
        pos_gene = pos.compress(to_keep)
        if len(G_gene) <= 1:
            # only a single mutation
            continue

        for dom_interval in dom[gene]:
            to_keep_dom = np.logical_and(
                pos_gene >= dom_interval[0], pos_gene <= dom_interval[1]
            )
            G_dom = G_gene.compress(to_keep_dom, axis=0)
            pos_dom = pos_gene.compress(to_keep_dom)
            if len(G_dom) <= 1:
                # only a single mutation in domain
                continue

            L, n = G_dom.shape
            G_dict, any_missing = moments.LD.Parsing._sparsify_genotype_matrix(G_dom)
            Counts = moments.LD.Parsing.spt.count_genotypes_sparse(
                G_dict, n, missing=any_missing
            )
            # each row of Counts is
            # [nAABB, nAABb, nAAbb, nAaBB, nAaBb, nAabb, naaBB, naaBb, naabb]
            genes_genotype_counts[gene].append(Counts)
            genes_positions[gene].append(pos_dom)
    return genes_genotype_counts, genes_positions


def get_two_locus_genotype_counts_outside_domains(genes, G, pos, dom):
    """
    For each gene, return genotype counts for each configuration seen.
    """
    genes_set = set(genes)
    genes_genotype_counts = defaultdict(list)
    genes_positions = defaultdict(list)
    for gene in genes_set:
        if gene not in dom.keys():
            continue

        to_keep = genes == gene
        G_gene = G.compress(to_keep, axis=0)
        pos_gene = pos.compress(to_keep)
        if len(G_gene) <= 1:
            # only a single mutation
            continue

        to_keep_outside_dom = np.ones(len(pos_gene), dtype=bool)
        for dom_interval in dom[gene]:
            to_keep_this_dom = np.logical_or(
                pos_gene < dom_interval[0], pos_gene > dom_interval[1]
            )
            to_keep_outside_dom = np.logical_and(to_keep_outside_dom, to_keep_this_dom)
        G_dom = G_gene.compress(to_keep_outside_dom, axis=0)
        pos_dom = pos_gene.compress(to_keep_outside_dom)
        if len(G_dom) <= 1:
            # only a single mutation in domain
            continue

        L, n = G_dom.shape
        G_dict, any_missing = moments.LD.Parsing._sparsify_genotype_matrix(G_dom)
        Counts = moments.LD.Parsing.spt.count_genotypes_sparse(
            G_dict, n, missing=any_missing
        )
        # each row of Counts is
        # [nAABB, nAABb, nAAbb, nAaBB, nAaBb, nAabb, naaBB, naaBb, naabb]
        genes_genotype_counts[gene].append(Counts)
        genes_positions[gene].append(pos_dom)
    return genes_genotype_counts, genes_positions


def compute_ld_stats(counts, nAmin=None, nAmax=None, nBmin=None, nBmax=None):
    nA = (
        2 * counts[:, 0]
        + 2 * counts[:, 1]
        + 2 * counts[:, 2]
        + counts[:, 3]
        + counts[:, 4]
        + counts[:, 5]
    )
    nB = (
        2 * counts[:, 0]
        + 2 * counts[:, 3]
        + 2 * counts[:, 6]
        + counts[:, 1]
        + counts[:, 4]
        + counts[:, 7]
    )
    keep = np.array([True for _ in nA])
    if nAmin is not None:
        keep[nA < nAmin] = False
    if nAmax is not None:
        keep[nA > nAmax] = False
    if nBmin is not None:
        keep[nB < nBmin] = False
    if nBmax is not None:
        keep[nB > nBmax] = False
    cs = counts.compress(keep, axis=0)
    if len(cs) == 0:
        return 0, [0, 0, 0]
    else:
        num_pairs = len(cs)
        D = moments.LD.Parsing.gcs.compute_D(cs)
        D2 = moments.LD.Parsing.gcs.compute_D2(cs)
        pi2 = moments.LD.Parsing.gcs.compute_pi2(cs)
        return num_pairs, [sum(D), sum(D2), sum(pi2)]


def compute_ld_stats_in_bins(
    all_counts, ds, b, nAmin=None, nAmax=None, nBmin=None, nBmax=None
):
    assert len(all_counts) == len(ds)
    to_keep = np.logical_and(ds > b[0], ds <= b[1])
    counts = all_counts.compress(to_keep, axis=0)
    if len(counts) == 0:
        return 0, [0, 0, 0]
    nA = (
        2 * counts[:, 0]
        + 2 * counts[:, 1]
        + 2 * counts[:, 2]
        + counts[:, 3]
        + counts[:, 4]
        + counts[:, 5]
    )
    nB = (
        2 * counts[:, 0]
        + 2 * counts[:, 3]
        + 2 * counts[:, 6]
        + counts[:, 1]
        + counts[:, 4]
        + counts[:, 7]
    )
    keep = np.array([True for _ in nA])
    if nAmin is not None:
        keep[nA < nAmin] = False
    if nAmax is not None:
        keep[nA > nAmax] = False
    if nBmin is not None:
        keep[nB < nBmin] = False
    if nBmax is not None:
        keep[nB > nBmax] = False
    counts = counts.compress(keep, axis=0)
    if len(counts) == 0:
        return 0, [0, 0, 0]
    else:
        num_pairs = len(counts)
        D = moments.LD.Parsing.gcs.compute_D(counts)
        D2 = moments.LD.Parsing.gcs.compute_D2(counts)
        pi2 = moments.LD.Parsing.gcs.compute_pi2(counts)
        return num_pairs, [sum(D), sum(D2), sum(pi2)]


def getOverlap(a, b):
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))


def bootstrap_genes(
    LD_by_gene, annots=["synonymous", "missense", "loss_of_function"], reps=1000
):
    SEs = {}
    for annot in annots:
        vals = np.array(LD_by_gene[annot])
        n = len(vals)
        replicates = []
        for i in range(reps):
            c = np.random.choice(n, size=n)
            new_vals = vals[c]
            D, D2, pi2 = new_vals.sum(axis=0)
            replicates.append([D / pi2, D2 / pi2])
        v = np.sqrt(np.var(replicates, axis=0))
        SEs[annot] = {"sd1": v[0], "sd2": v[1]}
    return SEs


def bootstrap_genes_bins(
    LD_by_gene_bins, annots=["synonymous", "missense", "loss_of_function"], reps=1000
):
    SEs = {}
    for annot in annots:
        SEs[annot] = {}
        vals = {
            b: np.array(LD_by_gene_bins[annot][b])
            for b in LD_by_gene_bins[annot].keys()
        }
        n = len(LD_by_gene_bins[annot][list(LD_by_gene_bins[annot].keys())[0]])
        replicates = {b: [] for b in LD_by_gene_bins[annot].keys()}
        for i in range(reps):
            c = np.random.choice(n, size=n)
            for b in LD_by_gene_bins[annot].keys():
                new_vals = vals[b][c]
                D, D2, pi2 = new_vals.sum(axis=0)
                replicates[b].append([D / pi2, D2 / pi2])
        for b in LD_by_gene_bins[annot].keys():
            v = np.sqrt(np.var(replicates[b], axis=0))
            SEs[annot][b] = {"sd1": v[0], "sd2": v[1]}
    return SEs


def compress_domains(ds):
    any_compressed = 1
    for d in ds:
        assert d[0] <= d[1]
    while any_compressed and len(ds) > 1:
        any_compressed = 0
        new_ds = [ds[0]]
        for d in ds[1:]:
            combined = 0
            for i, d2 in enumerate(new_ds):
                if getOverlap(d, d2) > 0:
                    new_ds[i] = (min(d[0], d2[0]), max(d[1], d2[1]))
                    any_compressed = 1
                    combined = 1
                    # combine and break
            if not combined:
                new_ds.append(d)
        ds = new_ds
    ds.sort()
    return ds


def merge_domains(domains):
    for chrom in domains.keys():
        for gene in domains[chrom].keys():
            ds = domains[chrom][gene]
            if len(ds) > 1:
                domains[chrom][gene] = compress_domains(ds)
    return domains


def pos_dists(pos):
    dists = np.zeros(len(pos) * (len(pos) - 1) // 2, dtype=int)
    c = 0
    for i, l in enumerate(pos[:-1]):
        for r in pos[i + 1 :]:
            dists[c] = r - l
            c += 1
    assert np.all(dists > 0)
    return dists


def get_distances(positions):
    distances = {}
    for chrom in positions.keys():
        distances[chrom] = {}
        for annot in positions[chrom].keys():
            distances[chrom][annot] = {}
            for gene in positions[chrom][annot].keys():
                distances[chrom][annot][gene] = []
                for pos in positions[chrom][annot][gene]:
                    distances[chrom][annot][gene].append(pos_dists(pos))
    return distances


def compute_LD_by_gene_list(
    genotype_counts, nAmin=None, nAmax=None, nBmin=None, nBmax=None
):
    """
    Given genotype counts (and any constraints on allele counts for loci A and B),
    return the number of pairs of mutations in total for each class of mutations,
    and the LD stats within each gene as a list, with each entry [D, D2, pi2],
    again for each class of mutations.
    """
    LD_by_gene = {annot: [] for annot in annots}
    num_pairs_by_gene = {annot: 0 for annot in annots}
    for chrom, annot_data in genotype_counts.items():
        for annot, gene_data in annot_data.items():
            for gene, count_list in gene_data.items():
                for counts in count_list:
                    num_pairs, gene_stats = compute_ld_stats(
                        counts, nAmin=nAmin, nAmax=nAmax, nBmin=nBmin, nBmax=nBmax
                    )
                    num_pairs_by_gene[annot] += num_pairs
                    LD_by_gene[annot].append(gene_stats)
    return num_pairs_by_gene, LD_by_gene


def compute_LD_by_gene_list_bins(
    genotype_counts,
    distances,
    bin_edges,
    nAmin=None,
    nAmax=None,
    nBmin=None,
    nBmax=None,
):
    LD_by_gene_bins = {annot: {b: [] for b in bin_edges} for annot in annots}
    num_pairs_by_gene_bins = {annot: {b: 0 for b in bin_edges} for annot in annots}
    for chrom, annot_data in genotype_counts.items():
        for annot, gene_data in annot_data.items():
            for gene, count_list in gene_data.items():
                ds = distances[chrom][annot][gene]
                assert len(count_list) == len(ds)
                for counts, d in zip(count_list, ds):
                    for b in bin_edges:
                        num_pairs, gene_stats = compute_ld_stats_in_bins(
                            counts, d, b, nAmin=nAmin, nAmax=nAmax, nBmin=nBmin, nBmax=nBmax
                        )
                        num_pairs_by_gene_bins[annot][b] += num_pairs
                        LD_by_gene_bins[annot][b].append(gene_stats)
    return num_pairs_by_gene_bins, LD_by_gene_bins


def get_matching_stats(
    LD_by_gene_outside_bins,
    num_pairs_by_gene_within_bins,
    bin_edges,
    annots=["synonymous", "missense"],
):
    data = {}
    for annot in annots:
        tot = sum(num_pairs_by_gene_within_bins[annot].values())
        weights = [num_pairs_by_gene_within_bins[annot][b] / tot for b in bin_edges]
        D, D2, pi2 = np.sum(
            [
                p * np.sum(LD_by_gene_outside_bins[annot][b], axis=0)
                for b, p in zip(bin_edges, weights)
            ],
            axis=0,
        )
        data[annot] = {"sd1": D / pi2, "sd2": D2 / pi2}
    return data


def bootstrap_matching(
    LD_by_gene_outside_bins,
    num_pairs_by_gene_within_bins,
    bin_edges,
    annots=["synonymous", "missense"],
    reps=1000,
):
    SEs = {}
    for annot in annots:
        n = len(
            LD_by_gene_outside_bins[annot][
                list(LD_by_gene_outside_bins[annot].keys())[0]
            ]
        )
        tot = sum(num_pairs_by_gene_within_bins[annot].values())
        weights = [num_pairs_by_gene_within_bins[annot][b] / tot for b in bin_edges]
        data = []
        for rep in range(reps):
            c = np.random.choice(n, size=n)
            D, D2, pi2 = np.sum(
                [
                    p * np.sum(np.array(LD_by_gene_outside_bins[annot][b])[c], axis=0)
                    for b, p in zip(bin_edges, weights)
                ],
                axis=0,
            )
            data.append([D / pi2, D2 / pi2])
        stderrs = np.sqrt(np.var(data, axis=0))
        SEs[annot] = {"sd1": stderrs[0], "sd2": stderrs[1]}
    return SEs


## save data
def get_data(LD_data):
    data_out = {}
    for annot in annots:
        data_out[annot] = {}
        D, D2, pi2 = np.sum(LD_data[annot], axis=0)
        data_out[annot]["sd1"] = D / pi2
        data_out[annot]["sd2"] = D2 / pi2
    return data_out


def get_data_binned(LD_data):
    data_out = {}
    for annot in annots:
        data_out[annot] = {"sd1": [], "sd2": []}
        for b in bin_edges:
            D, D2, pi2 = np.sum(LD_data[annot][b], axis=0)
            data_out[annot]["sd1"].append(D / pi2)
            data_out[annot]["sd2"].append(D2 / pi2)
    return data_out


if __name__ == "__main__":
    POP = sys.argv[1]

    domains = pickle.load(open("data/supp/domain_dict.bp", "rb"))
    domains = merge_domains(domains)

    annots = ["synonymous", "missense", "loss_of_function"]
    ld_win_dom = {ann: {} for ann in annots}
    num_win_dom = {ann: {} for ann in annots}

    ld_btw_dom = {ann: {} for ann in annots}
    num_btw_dom = {ann: {} for ann in annots}

    distances_within = {ann: [] for ann in annots}
    distances_between = {ann: [] for ann in annots}

    chroms = range(1, 23)
    annots = ["synonymous", "missense", "loss_of_function"]

    ## within domains
    genotype_counts = {}
    positions_within = {}
    for chrom in chroms:
        genotype_counts[chrom] = {}
        positions_within[chrom] = {}
        # print("processing chromosome", chrom)
        positions, annotations, genes, G = load_genotype_matrix(chrom, POP)
        # get stats within domains
        for annot in annots:
            pos_sub, genes_sub, G_sub = subset_annotation(
                positions, annotations, genes, G, annot=annot
            )
            cs, ps = get_two_locus_genotype_counts_within_domains(
                genes_sub, G_sub, pos_sub, domains[chrom]
            )
            genotype_counts[chrom][annot] = cs
            positions_within[chrom][annot] = ps

    # compute for all
    num_pairs_by_gene_within, LD_by_gene_within = compute_LD_by_gene_list(
        genotype_counts
    )

    # compute for n <= 2
    num_pairs_by_gene_within_leq_2, LD_by_gene_within_leq_2 = compute_LD_by_gene_list(
        genotype_counts, nAmax=2, nBmax=2
    )

    # compute for 3 <= n <= 8
    (
        num_pairs_by_gene_within_3_to_8,
        LD_by_gene_within_3_to_8,
    ) = compute_LD_by_gene_list(genotype_counts, nAmin=3, nAmax=8, nBmin=3, nBmax=8)

    # compute for n >= 9
    num_pairs_by_gene_within_geq_9, LD_by_gene_within_geq_9 = compute_LD_by_gene_list(
        genotype_counts, nAmin=9, nBmin=9
    )

    ## outside domains
    genotype_counts_outside = {}
    positions_outside = {}
    for chrom in chroms:
        genotype_counts_outside[chrom] = {}
        positions_outside[chrom] = {}
        # print("processing chromosome", chrom)
        positions, annotations, genes, G = load_genotype_matrix(chrom, POP)
        # get stats within domains
        for annot in annots:
            pos_sub, genes_sub, G_sub = subset_annotation(
                positions, annotations, genes, G, annot=annot
            )
            cs, ps = get_two_locus_genotype_counts_outside_domains(
                genes_sub, G_sub, pos_sub, domains[chrom]
            )
            genotype_counts_outside[chrom][annot] = cs
            positions_outside[chrom][annot] = ps

    # compute for all
    num_pairs_by_gene_outside, LD_by_gene_outside = compute_LD_by_gene_list(
        genotype_counts_outside
    )

    # compute for n <= 2
    num_pairs_by_gene_outside_leq_2, LD_by_gene_outside_leq_2 = compute_LD_by_gene_list(
        genotype_counts_outside, nAmax=2, nBmax=2
    )

    # compute for 3 <= n <= 8
    (
        num_pairs_by_gene_outside_3_to_8,
        LD_by_gene_outside_3_to_8,
    ) = compute_LD_by_gene_list(
        genotype_counts_outside, nAmin=3, nAmax=8, nBmin=3, nBmax=8
    )

    # compute for n >= 9
    num_pairs_by_gene_outside_geq_9, LD_by_gene_outside_geq_9 = compute_LD_by_gene_list(
        genotype_counts_outside, nAmin=9, nBmin=9
    )

    ## outside domains, matched to within domain distances
    distances_within = get_distances(positions_within)
    distances_outside = get_distances(positions_outside)

    bins = np.array(
        [0, 100, 200, 300, 500, 1000, 2000, 3000, 5000, 10000, 20000, 50000, 5e6]
    )
    bin_edges = list(zip(bins[:-1], bins[1:]))

    # all pairs within domains by distance
    (
        num_pairs_by_gene_within_bins,
        LD_by_gene_within_bins,
    ) = compute_LD_by_gene_list_bins(
        genotype_counts,
        distances_within,
        bin_edges,
    )

    # pairs n <= 2 within domains by distances
    (
        num_pairs_by_gene_within_bins_leq_2,
        LD_by_gene_within_bins_leq_2,
    ) = compute_LD_by_gene_list_bins(
        genotype_counts,
        distances_within,
        bin_edges,
        nAmax=2,
        nBmax=2,
    )

    # pairs 3 <= n <= 8 within domains by distances
    (
        num_pairs_by_gene_within_bins_3_to_8,
        LD_by_gene_within_bins_3_to_8,
    ) = compute_LD_by_gene_list_bins(
        genotype_counts,
        distances_within,
        bin_edges,
        nAmin=3,
        nAmax=8,
        nBmin=3,
        nBmax=8,
    )

    # pairs n >= 9 within domains by distances
    (
        num_pairs_by_gene_within_bins_geq_9,
        LD_by_gene_within_bins_geq_9,
    ) = compute_LD_by_gene_list_bins(
        genotype_counts,
        distances_within,
        bin_edges,
        nAmin=9,
        nBmin=9,
    )

    # all pairs outside domains by distance
    (
        num_pairs_by_gene_outside_bins,
        LD_by_gene_outside_bins,
    ) = compute_LD_by_gene_list_bins(
        genotype_counts_outside,
        distances_outside,
        bin_edges,
    )

    # pairs n <= 2 within domains by distances
    (
        num_pairs_by_gene_outside_bins_leq_2,
        LD_by_gene_outside_bins_leq_2,
    ) = compute_LD_by_gene_list_bins(
        genotype_counts_outside,
        distances_outside,
        bin_edges,
        nAmax=2,
        nBmax=2,
    )

    # pairs 3 <= n <= 8 within domains by distances
    (
        num_pairs_by_gene_outside_bins_3_to_8,
        LD_by_gene_outside_bins_3_to_8,
    ) = compute_LD_by_gene_list_bins(
        genotype_counts_outside,
        distances_outside,
        bin_edges,
        nAmin=3,
        nAmax=8,
        nBmin=3,
        nBmax=8,
    )

    # pairs n >= 9 within domains by distances
    (
        num_pairs_by_gene_outside_bins_geq_9,
        LD_by_gene_outside_bins_geq_9,
    ) = compute_LD_by_gene_list_bins(
        genotype_counts_outside,
        distances_outside,
        bin_edges,
        nAmin=9,
        nBmin=9,
    )

    ## get all domain (within and outside) and save
    data = {}
    data["within"] = {
        "all": get_data(LD_by_gene_within),
        "leq2": get_data(LD_by_gene_within_leq_2),
        "3to8": get_data(LD_by_gene_within_3_to_8),
        "geq9": get_data(LD_by_gene_within_geq_9),
    }
    data["outside"] = {}
    data["outside"]["all_pairs"] = {
        "all": get_data(LD_by_gene_outside),
        "leq2": get_data(LD_by_gene_outside_leq_2),
        "3to8": get_data(LD_by_gene_outside_3_to_8),
        "geq9": get_data(LD_by_gene_outside_geq_9),
    }
    data["outside"]["matching"] = {
        "all": get_matching_stats(
            LD_by_gene_outside_bins, num_pairs_by_gene_within_bins, bin_edges
        ),
        "leq2": get_matching_stats(
            LD_by_gene_outside_bins_leq_2,
            num_pairs_by_gene_within_bins_leq_2,
            bin_edges,
        ),
        "3to8": get_matching_stats(
            LD_by_gene_outside_bins_3_to_8,
            num_pairs_by_gene_within_bins_3_to_8,
            bin_edges,
        ),
        "geq9": get_matching_stats(
            LD_by_gene_outside_bins_geq_9,
            num_pairs_by_gene_within_bins_geq_9,
            bin_edges,
        ),
    }
    data["bins"] = {
        "within": get_data_binned(LD_by_gene_within_bins),
        "outside": get_data_binned(LD_by_gene_outside_bins),
    }
    data["bin_edges"] = bin_edges

    pickle.dump(data, open(f"parsed_data_v2/{POP}.unphased.domains.bp", "wb+"))

    # get standard errors
    SEs = {}
    SEs["within"] = {
        "all": bootstrap_genes(LD_by_gene_within),
        "leq2": bootstrap_genes(LD_by_gene_within_leq_2),
        "3to8": bootstrap_genes(LD_by_gene_within_3_to_8),
        "geq9": bootstrap_genes(LD_by_gene_within_geq_9),
    }
    SEs["outside"] = {}
    SEs["outside"]["all_pairs"] = {
        "all": bootstrap_genes(LD_by_gene_outside),
        "leq2": bootstrap_genes(LD_by_gene_outside_leq_2),
        "3to8": bootstrap_genes(LD_by_gene_outside_3_to_8),
        "geq9": bootstrap_genes(LD_by_gene_outside_geq_9),
    }
    SEs["outside"]["matching"] = {
        "all": bootstrap_matching(
            LD_by_gene_outside_bins, num_pairs_by_gene_within_bins, bin_edges
        ),
        "leq2": bootstrap_matching(
            LD_by_gene_outside_bins_leq_2,
            num_pairs_by_gene_within_bins_leq_2,
            bin_edges,
        ),
        "3to8": bootstrap_matching(
            LD_by_gene_outside_bins_3_to_8,
            num_pairs_by_gene_within_bins_3_to_8,
            bin_edges,
        ),
        "geq9": bootstrap_matching(
            LD_by_gene_outside_bins_geq_9,
            num_pairs_by_gene_within_bins_geq_9,
            bin_edges,
        ),
    }

    pickle.dump(SEs, open(f"parsed_data_v2/{POP}.unphased.domains.SEs.bp", "wb+"))
