import numpy as np
import gzip
import moments.LD
import sys
import pickle


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


def compute_LD_within_genes(genes, G):
    """
    Compute LD statistics within each gene for the given genotype matrix.
    Also returns the number of pair comparisons made per gene.
    """
    genes_set = set(genes)
    gene_stat_sums = {}
    num_pairs_per_gene = {}
    for gene in genes_set:
        to_keep = genes == gene
        G_gene = G.compress(to_keep, axis=0)
        if len(G_gene) > 1:
            D2, Dz, pi2, D = moments.LD.Parsing.compute_pairwise_stats(
                G_gene, genotypes=True
            )
            gene_stat_sums[gene] = np.array(
                [np.sum(D2), np.sum(Dz), np.sum(pi2), np.sum(D)]
            )
            num_pairs_per_gene[gene] = len(D2)
        else:
            gene_stat_sums[gene] = np.zeros(4)
            num_pairs_per_gene[gene] = 0
    return gene_stat_sums, num_pairs_per_gene


def compute_LD_within_genes_by_frequency(genes, G, n_min=2, n_max=2):
    genes_set = set(genes)
    gene_stat_sums = {}
    num_pairs_per_gene = {}
    allele_counts = G.sum(axis=1)
    to_keep_ac = np.logical_and(allele_counts >= n_min, allele_counts <= n_max)
    for gene in genes_set:
        to_keep = genes == gene
        to_keep = np.logical_and(to_keep, to_keep_ac)
        G_gene = G.compress(to_keep, axis=0)
        if len(G_gene) > 1:
            D2, Dz, pi2, D = moments.LD.Parsing.compute_pairwise_stats(
                G_gene, genotypes=True
            )
            gene_stat_sums[gene] = np.array(
                [np.sum(D2), np.sum(Dz), np.sum(pi2), np.sum(D)]
            )
            num_pairs_per_gene[gene] = len(D2)
        else:
            gene_stat_sums[gene] = np.zeros(4)
            num_pairs_per_gene[gene] = 0
    return gene_stat_sums, num_pairs_per_gene


def compute_LD_by_distance(
    positions, genes, G, bp_bins=[0, 1e4, 2e4, 3e4, 4e4, 5e4], exclude_within=True
):
    """
    Compute LD only for pairs lying in different genes, within the given bp window
    """
    genes_set = set(genes)
    bins = [(l, r) for l, r in zip(bp_bins[:-1], bp_bins[1:])]
    stat_sums = {b: np.zeros(4) for b in bins}
    num_pairs = {b: 0 for b in bins}
    for ii, (pos, gene) in enumerate(zip(positions, genes)):
        other_gene = genes != gene
        for b in bins:
            to_keep = np.logical_and(
                positions - positions[ii] > b[0], positions - positions[ii] <= b[1]
            )
            if exclude_within:
                to_keep = np.logical_and(to_keep, other_gene)
            focal_haplotype = G[ii]
            compare_haplotypes = G.compress(to_keep, axis=0)
            if len(compare_haplotypes) > 0:
                D2, Dz, pi2, D = moments.LD.Parsing.compute_pairwise_stats_between(
                    focal_haplotype[np.newaxis, :], compare_haplotypes, genotypes=True
                )
                stat_sums[b] += np.array(
                    [np.sum(D2), np.sum(Dz), np.sum(pi2), np.sum(D)]
                )
                num_pairs[b] += len(compare_haplotypes)
    return stat_sums, num_pairs


if __name__ == "__main__":
    POP = sys.argv[1]

    ld_stats = {"synonymous": {}, "missense": {}, "loss_of_function": {}}
    num_pairs = {"synonymous": {}, "missense": {}, "loss_of_function": {}}

    bins = [0, 50e3, 100e3, 150e3, 200e3, 300e3, 500e3]
    ld_stats_btw = {
        "synonymous": {(l, r): np.zeros(4) for l, r in zip(bins[:-1], bins[1:])},
        "missense": {(l, r): np.zeros(4) for l, r in zip(bins[:-1], bins[1:])},
        "loss_of_function": {(l, r): np.zeros(4) for l, r in zip(bins[:-1], bins[1:])},
    }
    num_pairs_btw = {
        "synonymous": {(l, r): 0 for l, r in zip(bins[:-1], bins[1:])},
        "missense": {(l, r): 0 for l, r in zip(bins[:-1], bins[1:])},
        "loss_of_function": {(l, r): 0 for l, r in zip(bins[:-1], bins[1:])},
    }

    chroms = range(1, 23)

    for chrom in chroms:
        print("processing chromosome", chrom)
        positions, annotations, genes, G = load_genotype_matrix(chrom, POP)
        for annot in ld_stats.keys():
            pos_sub, genes_sub, G_sub = subset_annotation(
                positions, annotations, genes, G, annot=annot
            )
            LD, N = compute_LD_within_genes(genes_sub, G_sub)
            ld_stats[annot].update(LD)
            num_pairs[annot].update(N)
            LD_btw, N_btw = compute_LD_by_distance(
                pos_sub, genes_sub, G_sub, bp_bins=bins
            )
            for k, v in LD_btw.items():
                ld_stats_btw[annot][k] += v
            for k, v in N_btw.items():
                num_pairs_btw[annot][k] += v

    ld_stats_within = {
        annot: np.sum(list(ld_stats[annot].values()), axis=0)
        for annot in ld_stats.keys()
    }
    num_pairs_within = {
        annot: sum(num_pairs[annot].values()) for annot in num_pairs.keys()
    }

    outputs = {
        "ld_within": ld_stats_within,
        "num_pairs_within": num_pairs_within,
        "ld_between": ld_stats_btw,
        "num_pairs_between": num_pairs_btw,
        "bins": bins,
    }


    pickle.dump(outputs, open(f"parsed_data/{POP}.unphased.within.between.bp", "wb+"))

    """
    ld_stats_by_freqs = {}
    num_pairs_by_freqs = {}
    for m, M in [(1, 1), (2, 2), (3, 4), (5, 8), (9, 16), (17, 32), (64, 108)]:
        print(m, M)
        ld_stats_by_freqs[(m, M)] = {"synonymous": {}, "missense": {}, "loss_of_function": {}}
        num_pairs_by_freqs[(m, M)] = {"synonymous": {}, "missense": {}, "loss_of_function": {}}
        chroms = range(1, 23)
        for chrom in chroms:
            positions, annotations, genes, G = load_genotype_matrix(chrom, POP)
            for annot in ld_stats.keys():
                pos_sub, genes_sub, G_sub = subset_annotation(
                    positions, annotations, genes, G, annot=annot)
                LD, N = compute_LD_within_genes_by_frequency(genes_sub, G_sub, n_min=m, n_max=M)
                ld_stats_by_freqs[(m, M)][annot].update(LD)
                num_pairs_by_freqs[(m, M)][annot].update(N)
    """
