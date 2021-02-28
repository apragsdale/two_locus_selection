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


def compute_LD_outside_domains(pos, genes, G, dom, bins):
    """
    Compute LD statistics within each gene for the given genotype matrix.
    Also returns the number of pair comparisons made per gene.
    Only keep pairs of SNPs that fall within the same domain in each gene.
    """
    genes_set = set(genes)
    gene_stat_sums = {}
    num_pairs_per_gene = {}
    for gene in genes_set:
        to_keep = genes == gene
        G_gene = G.compress(to_keep, axis=0)
        pos_gene = pos.compress(to_keep)
        gene_stat_sums[gene] = {}
        num_pairs_per_gene[gene] = {}
        added_gene = False
        if len(G_gene) > 1:
            to_keep_dom = np.array([True for _ in pos_gene])
            for dom_interval in dom[gene]:
                to_keep_dom = np.logical_and(
                    to_keep_dom,
                    np.logical_or(
                        pos_gene < dom_interval[0], pos_gene > dom_interval[1]
                    ),
                )
            G_dom = G_gene.compress(to_keep_dom, axis=0)
            pos_dom = pos_gene.compress(to_keep_dom)
            if len(G_dom) > 1:
                D2, Dz, pi2, D = moments.LD.Parsing.compute_pairwise_stats(
                    G_dom, genotypes=True
                )
                # get pair distances
                pair_distances = []
                for i, l in enumerate(pos_dom[:-1]):
                    for r in pos_dom[i + 1 :]:
                        pair_distances.append(r - l)
                pair_distances = np.array(pair_distances)
                assert np.all(pair_distances > 0)
                # sum statistics over bins
                for b in bins:
                    added_gene = True
                    to_keep_bin = np.logical_and(
                        pair_distances >= b[0], pair_distances < b[1]
                    )
                    D2_bin = D2.compress(to_keep_bin)
                    Dz_bin = Dz.compress(to_keep_bin)
                    pi2_bin = pi2.compress(to_keep_bin)
                    D_bin = D.compress(to_keep_bin)
                    gene_stat_sums[gene][b] = np.array(
                        [np.sum(D2_bin), np.sum(Dz_bin), np.sum(pi2_bin), np.sum(D_bin)]
                    )
                    num_pairs_per_gene[gene][b] = np.sum(to_keep_bin)
        if not added_gene:
            gene_stat_sums[gene] = {b: np.zeros(4) for b in bins}
            num_pairs_per_gene[gene] = {b: 0 for b in bins}
    return gene_stat_sums, num_pairs_per_gene


def sum_over_genes(values, bins):
    data = {b: np.zeros(4) for b in bins}
    for v in values:
        for k, d in v.items():
            data[k] += d
    return data


def sum_over_weights(ld_stats, num_pairs, bin_weights):
    weighted_stats = np.zeros(4)
    for ii, b in enumerate(ld_stats.keys()):
        if num_pairs[b] > 0:
            weighted_stats += ld_stats[b] / num_pairs[b] * bin_weights[ii]
    return weighted_stats


def bootstrap_over_genes(
    gene_data, num_pairs, gene_sets, bins, bin_weights=None, reps=1000
):
    bs_sums = {k: {b: np.zeros(4) for b in bins} for k in gene_sets["gene_sets"].keys()}
    bs_nums = {k: {b: 0 for b in bins} for k in gene_sets["gene_sets"].keys()}
    for gene in gene_data:
        try:
            k = gene_sets["gene_set_map"][gene]
        except KeyError:
            # not good...
            k = 500
        for b in bins:
            bs_sums[k][b] += gene_data[gene][b]
            bs_nums[k][b] += num_pairs[gene][b]

    sigma_d1s = []
    for i in range(reps):
        choices = np.random.choice(range(len(bs_sums)), len(bs_sums))
        ld_stats_rep = {b: np.zeros(4) for b in bins}
        num_pairs_rep = {b: 0 for b in bins}
        for c in choices:
            for b in bins:
                ld_stats_rep[b] += bs_sums[c][b]
                num_pairs_rep[b] += bs_nums[c][b]
        if bin_weights is not None:
            bs_stats = sum_over_weights(ld_stats_rep, num_pairs_rep, bin_weights)
        else:
            bs_stats = np.sum(list(ld_stats_rep.values()), axis=0)
        sigma_d1s.append(bs_stats[3] / bs_stats[2])
    return np.std(sigma_d1s)


if __name__ == "__main__":
    POP = sys.argv[1]

    domains = pickle.load(open("data/supp/domain_dict.bp", "rb"))

    # average domain distances
    distance_data = pickle.load(
        open(f"parsed_data/{POP}.distances_within_between_domains.bp", "rb")
    )
    bin_edges = distance_data["bin_edges"]
    bins = [(l, r) for l, r in zip(bin_edges[:-1], bin_edges[1:])]
    num_within_domains = distance_data["num_within"]
    num_between_domains = distance_data["num_between"]

    annots = ["synonymous", "missense", "loss_of_function"]
    ld_outside = {ann: {} for ann in annots}
    num_outside = {ann: {} for ann in annots}

    chroms = range(1, 23)
    for chrom in chroms:
        #print("processing chromosome", chrom)
        positions, annotations, genes, G = load_genotype_matrix(chrom, POP)
        # get stats within domains
        for annot in ld_outside.keys():
            pos_sub, genes_sub, G_sub = subset_annotation(
                positions, annotations, genes, G, annot=annot
            )
            LD, N = compute_LD_outside_domains(
                pos_sub, genes_sub, G_sub, domains[chrom], bins
            )
            ld_outside[annot].update(LD)
            num_outside[annot].update(N)

    ld_stats = {
        annot: sum_over_genes(ld_outside[annot].values(), bins)
        for annot in ld_outside.keys()
    }
    num_pairs = {
        annot: {
            b: sum([num_outside[annot][g][b] for g in num_outside[annot].keys()])
            for b in bins
        }
        for annot in num_outside.keys()
    }

    ld_outside_totals = {
        annot: np.sum(list(ld_stats[annot].values()), axis=0) for annot in annots
    }
    ld_weighted_within = {
        annot: sum_over_weights(
            ld_stats[annot], num_pairs[annot], num_within_domains[annot][0]
        )
        for annot in annots
    }
    ld_weighted_between = {
        annot: sum_over_weights(
            ld_stats[annot], num_pairs[annot], num_between_domains[annot][0]
        )
        for annot in annots
    }

    outputs = {
        "ld_outside_totals": ld_outside_totals,
        "ld_weighted_within": ld_weighted_within,
        "ld_weighted_between": ld_weighted_between,
    }

    pickle.dump(outputs, open(f"parsed_data/{POP}.unphased.outside_domains.bp", "wb+"))

    print(POP)
    # bootstrap sigma_d^1 stats
    gene_sets = pickle.load(open("data/supp/bootstrap_gene_sets.500.bp", "rb"))
    # outside totals:
    bs_outside = {}
    for annot in ld_outside.keys():
        bs_outside[annot] = bootstrap_over_genes(
            ld_outside[annot], num_outside[annot], gene_sets, bins
        )
        print(
            annot,
            ld_outside_totals[annot][3] / ld_outside_totals[annot][2],
            bs_outside[annot],
        )

    bs_within = {}
    for annot in ld_outside.keys():
        bs_within[annot] = bootstrap_over_genes(
            ld_outside[annot],
            num_outside[annot],
            gene_sets,
            bins,
            bin_weights=num_within_domains[annot][0],
        )
        print(
            annot,
            ld_weighted_within[annot][3] / ld_weighted_within[annot][2],
            bs_within[annot],
        )

    bs_between = {}
    for annot in ld_outside.keys():
        bs_between[annot] = bootstrap_over_genes(
            ld_outside[annot],
            num_outside[annot],
            gene_sets,
            bins,
            bin_weights=num_between_domains[annot][0],
        )
        print(
            annot,
            ld_weighted_between[annot][3] / ld_weighted_between[annot][2],
            bs_between[annot],
        )

    bs_outputs = {
        "bs_outside_domains": bs_outside,
        "bs_weighted_within": bs_within,
        "bs_weighted_between": bs_between,
    }

    pickle.dump(
        bs_outputs,
        open(
            f"parsed_data/{POP}.unphased._outside_domains.sigmad1.bootstrap_std_err.bp",
            "wb+",
        ),
    )
