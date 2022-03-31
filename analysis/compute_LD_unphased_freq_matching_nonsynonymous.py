import numpy as np
import gzip
import moments.LD
import sys
import pickle

"""
Combine missense and lof into "nonsynonymous"
"""

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
    if annot not in ["synonymous", "nonsynonymous"]:
        raise ValueError("pick one of synonymous, nonsynonymous")
    to_keep = np.zeros(len(positions))
    for ii, pos_annot in enumerate(annotations):
        is_non = False
        is_syn = False
        for lof_annot in lof_annotations:
            if pos_annot.startswith(lof_annot):
                is_non = True
        if not is_non:
            if pos_annot.startswith("missense"):
                is_non = True
        if not is_non:
            if pos_annot.startswith("synonymous"):
                is_syn = True
        if (annot == "nonsynonymous") and is_non:
            to_keep[ii] = 1
        elif annot == "synonymous" and is_syn:
            to_keep[ii] = 1

    new_positions = positions.compress(to_keep)
    new_genes = genes.compress(to_keep)
    new_G = G.compress(to_keep, axis=0)
    return new_positions, new_genes, new_G


def get_two_locus_genotype_counts(genes, G, pos):
    """
    For each gene, return genotype counts for each configuration seen.
    """
    genes_set = set(genes)
    genes_genotype_counts = {}
    genes_positions = {}
    for gene in genes_set:
        to_keep = genes == gene
        G_gene = G.compress(to_keep, axis=0)
        pos_gene = pos.compress(to_keep)
        if len(G_gene) <= 1:
            # only a single mutation
            continue
        L, n = G_gene.shape
        G_dict, any_missing = moments.LD.Parsing._sparsify_genotype_matrix(G_gene)
        Counts = moments.LD.Parsing.spt.count_genotypes_sparse(
            G_dict, n, missing=any_missing
        )
        # each row of Counts is
        # [nAABB, nAABb, nAAbb, nAaBB, nAaBb, nAabb, naaBB, naaBb, naabb]
        genes_genotype_counts[gene] = Counts
        genes_positions[gene] = pos_gene
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


def bootstrap_genes(
    LD_by_gene, annots=["synonymous", "nonsynonymous"], reps=1000
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
    LD_by_gene_bins, annots=["synonymous", "nonsynonymous"], reps=1000
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
                distances[chrom][annot][gene] = pos_dists(positions[chrom][annot][gene])
    return distances


def compute_LD_by_gene(genotype_counts, nAmin=None, nAmax=None, nBmin=None, nBmax=None):
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
            for gene, counts in gene_data.items():
                num_pairs, gene_stats = compute_ld_stats(
                    counts, nAmin=nAmin, nAmax=nAmax, nBmin=nBmin, nBmax=nBmax
                )
                num_pairs_by_gene[annot] += num_pairs
                LD_by_gene[annot].append(gene_stats)
    return num_pairs_by_gene, LD_by_gene

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

    chroms = range(1, 23)
    annots = ["synonymous", "nonsynonymous"]
    # get genotype counts for each gene in all chromosomes
    genotype_counts = {}
    positions_in_genes = {}
    for chrom in chroms:
        genotype_counts[chrom] = {}
        positions_in_genes[chrom] = {}
        # print("processing chromosome", chrom)
        positions, annotations, genes, G = load_genotype_matrix(chrom, POP)
        for annot in annots:
            pos_sub, genes_sub, G_sub = subset_annotation(
                positions, annotations, genes, G, annot=annot
            )
            cs, ps = get_two_locus_genotype_counts(
                genes_sub,
                G_sub,
                pos_sub,
            )
            genotype_counts[chrom][annot] = cs
            positions_in_genes[chrom][annot] = ps

    # compute for all
    num_pairs_by_gene, LD_by_gene = compute_LD_by_gene(genotype_counts)

    # compute for n <= 2
    num_pairs_by_gene_leq_2, LD_by_gene_leq_2 = compute_LD_by_gene(
        genotype_counts, nAmax=2, nBmax=2
    )

    # compute for 3 <= n <= 8
    num_pairs_by_gene_3_to_8, LD_by_gene_3_to_8 = compute_LD_by_gene(
        genotype_counts, nAmin=3, nBmin=3, nAmax=8, nBmax=8
    )

    # compute for n >= 9
    num_pairs_by_gene_geq_9, LD_by_gene_geq_9 = compute_LD_by_gene(
        genotype_counts, nAmin=9, nBmin=9
    )

    # by bins
    distances = get_distances(positions_in_genes)

    bins = np.array(
        [0, 100, 200, 300, 500, 1000, 2000, 3000, 5000, 10000, 20000, 50000, 5e6]
    )
    bin_edges = list(zip(bins[:-1], bins[1:]))

    LD_by_gene_bins = {annot: {b: [] for b in bin_edges} for annot in annots}
    num_pairs_by_gene_bins = {annot: {b: 0 for b in bin_edges} for annot in annots}
    for chrom, annot_data in genotype_counts.items():
        for annot, gene_data in annot_data.items():
            for gene, counts in gene_data.items():
                ds = distances[chrom][annot][gene]
                assert len(ds) == len(counts)
                for b in bin_edges:
                    num_pairs, gene_stats = compute_ld_stats_in_bins(
                        counts, ds, b, nAmin=None, nAmax=None, nBmin=None, nBmax=None
                    )
                    num_pairs_by_gene_bins[annot][b] += num_pairs
                    LD_by_gene_bins[annot][b].append(gene_stats)

    ## save data
    data = {}
    data["all_sites"] = {
        "all": get_data(LD_by_gene),
        "leq2": get_data(LD_by_gene_leq_2),
        "3to8": get_data(LD_by_gene_3_to_8),
        "geq9": get_data(LD_by_gene_geq_9),
    }
    data["bins"] = {
        "all_sites": get_data_binned(LD_by_gene_bins),
    }
    data["bin_edges"] = bin_edges

    pickle.dump(data, open(f"parsed_data_v2/{POP}.unphased.nonsyn.all.bp", "wb+"))

    SEs = {}
    SEs["all_sites"] = {
        "all": bootstrap_genes(LD_by_gene),
        "leq2": bootstrap_genes(LD_by_gene_leq_2),
        "3to8": bootstrap_genes(LD_by_gene_3_to_8),
        "geq9": bootstrap_genes(LD_by_gene_geq_9),
    }
    SEs["bins"] = {"all_sites": bootstrap_genes_bins(LD_by_gene_bins)}

    pickle.dump(SEs, open(f"parsed_data_v2/{POP}.unphased.nonsyn.all.SEs.bp", "wb+"))
