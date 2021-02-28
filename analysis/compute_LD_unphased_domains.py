import numpy as np
import itertools
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


def compute_LD_within_domains(pos, genes, G, dom):
    """
    Compute LD statistics within each gene for the given genotype matrix.
    Also returns the number of pair comparisons made per gene.
    Only keep pairs of SNPs that fall within the same domain in each gene.
    """
    genes_set = set(genes)
    gene_stat_sums = {}
    num_pairs_per_gene = {}
    found = 0
    not_found = 0
    distances = []
    for gene in genes_set:
        if gene not in dom.keys():
            # print(gene, "not found")
            not_found += 1
            continue
        else:
            found += 1
        to_keep = genes == gene
        G_gene = G.compress(to_keep, axis=0)
        pos_gene = pos.compress(to_keep)
        gene_stat_sums[gene] = np.zeros(4)
        num_pairs_per_gene[gene] = 0
        if len(G_gene) > 1:
            for dom_interval in dom[gene]:
                to_keep_dom = np.logical_and(
                    pos_gene >= dom_interval[0], pos_gene <= dom_interval[1]
                )
                G_dom = G_gene.compress(to_keep_dom, axis=0)
                if len(G_dom) > 1:
                    D2, Dz, pi2, D = moments.LD.Parsing.compute_pairwise_stats(
                        G_dom, genotypes=True
                    )
                    gene_stat_sums[gene] += np.array(
                        [np.sum(D2), np.sum(Dz), np.sum(pi2), np.sum(D)]
                    )
                    num_pairs_per_gene[gene] += len(D2)
                    pos_dom = pos_gene.compress(to_keep_dom)
                    distances_dom = [
                        abs(r - l) for l, r in itertools.combinations(pos_dom, 2)
                    ]
                    distances += distances_dom
    print("found:", found, ", not found:", not_found)
    return gene_stat_sums, num_pairs_per_gene, distances


def getOverlap(a, b):
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))


def compute_LD_between_domains(positions, genes, G, dom):
    """
    Compute LD statistics within each gene for the given genotype matrix.
    Also returns the number of pair comparisons made per gene.
    Only keep pairs of SNPs that fall within a domain but different domains.
    """
    genes_set = set(genes)
    gene_stat_sums = {}
    num_pairs_per_gene = {}
    distances = []
    for ii, (pos, gene) in enumerate(zip(positions, genes)):
        if gene not in dom.keys():
            continue
        to_keep = genes == gene
        G_gene = G.compress(to_keep, axis=0)
        pos_gene = positions.compress(to_keep)
        gene_stat_sums.setdefault(gene, np.zeros(4))
        num_pairs_per_gene.setdefault(gene, 0)
        if len(G_gene) > 1:
            in_domain = False
            for dom_interval in dom[gene]:
                if pos >= dom_interval[0] and pos <= dom_interval[1]:
                    pos_domain = dom_interval
                    in_domain = True
                    break
            if in_domain:
                focal_G = G[ii]
                to_keep_dom = np.zeros(len(pos_gene))
                # only keep other domains that do not overlap
                other_domains = []
                for other_dom_interval in dom[gene]:
                    if getOverlap(pos_domain, other_dom_interval) <= 0:
                        other_domains.append(other_dom_interval)
                if len(other_domains) > 0:
                    for jj, other_pos in enumerate(pos_gene):
                        if other_pos <= pos:
                            continue
                        else:
                            for other_domain in other_domains:
                                if (
                                    other_pos >= other_domain[0]
                                    and other_pos <= other_domain[1]
                                ):
                                    to_keep_dom[jj] = 1
                    # keep positions greater than focal pos and in other domain
                    G_dom = G_gene.compress(to_keep_dom, axis=0)
                    pos_dom = pos_gene.compress(to_keep_dom)
                    if len(G_dom) >= 1:
                        (
                            D2,
                            Dz,
                            pi2,
                            D,
                        ) = moments.LD.Parsing.compute_pairwise_stats_between(
                            focal_G[np.newaxis, :], G_dom, genotypes=True
                        )
                        gene_stat_sums[gene] += np.array(
                            [np.sum(D2), np.sum(Dz), np.sum(pi2), np.sum(D)]
                        )
                        num_pairs_per_gene[gene] += len(D2)
                        distances += [abs(pos - other_pos) for other_pos in pos_dom]
    return gene_stat_sums, num_pairs_per_gene, distances


def bootstrap_over_genes(gene_data, gene_sets, reps=1000):
    bs_sums = [np.zeros(4) for k in gene_sets["gene_sets"].keys()]
    for gene in gene_data:
        try:
            bs_sums[gene_sets["gene_set_map"][gene]] += gene_data[gene]
        except KeyError:
            print(gene)
            # not good...
            bs_sums[-1] += gene_data[gene]

    sigma_d1s = []
    for i in range(reps):
        choices = np.random.choice(range(len(bs_sums)), len(bs_sums))
        bs_stats = np.sum([bs_sums[j] for j in choices], axis=0)
        sigma_d1s.append(bs_stats[3] / bs_stats[2])
    return np.std(sigma_d1s)


def bootstrap_over_genes_meanD(gene_data, num_pairs, gene_sets, reps=1000):
    bs_sums = [np.zeros(4) for k in gene_sets["gene_sets"].keys()]
    num_sums = [0 for k in gene_sets["gene_sets"].keys()]
    for gene in gene_data:
        try:
            bs_sums[gene_sets["gene_set_map"][gene]] += gene_data[gene]
            num_sums[gene_sets["gene_set_map"][gene]] += num_pairs[gene]
        except KeyError:
            print(gene)
            # not good...
            bs_sums[-1] += gene_data[gene]
            num_sums[-1] += num_pairs[gene]

    mean_Ds = []
    for i in range(reps):
        choices = np.random.choice(range(len(bs_sums)), len(bs_sums))
        bs_stats = np.sum([bs_sums[j] for j in choices], axis=0)
        total_pairs = np.sum([num_sums[j] for j in choices])
        mean_Ds.append(bs_stats[3] / total_pairs)
    return np.std(mean_Ds)


if __name__ == "__main__":
    POP = sys.argv[1]

    domains = pickle.load(open("data/supp/domain_dict.bp", "rb"))

    annots = ["synonymous", "missense", "loss_of_function"]
    ld_win_dom = {ann: {} for ann in annots}
    num_win_dom = {ann: {} for ann in annots}

    ld_btw_dom = {ann: {} for ann in annots}
    num_btw_dom = {ann: {} for ann in annots}

    distances_within = {ann: [] for ann in annots}
    distances_between = {ann: [] for ann in annots}

    chroms = range(1, 23)
    for chrom in chroms:
        print("processing chromosome", chrom)
        positions, annotations, genes, G = load_genotype_matrix(chrom, POP)
        # get stats within domains
        for annot in ld_win_dom.keys():
            pos_sub, genes_sub, G_sub = subset_annotation(
                positions, annotations, genes, G, annot=annot
            )
            LD, N, d = compute_LD_within_domains(
                pos_sub, genes_sub, G_sub, domains[chrom]
            )
            ld_win_dom[annot].update(LD)
            num_win_dom[annot].update(N)
            distances_within[annot] += d
        for annot in ld_btw_dom.keys():
            pos_sub, genes_sub, G_sub = subset_annotation(
                positions, annotations, genes, G, annot=annot
            )
            LD, N, d = compute_LD_between_domains(
                pos_sub, genes_sub, G_sub, domains[chrom]
            )
            ld_btw_dom[annot].update(LD)
            num_btw_dom[annot].update(N)
            distances_between[annot] += d

    ld_stats_within = {
        annot: np.sum(list(ld_win_dom[annot].values()), axis=0)
        for annot in ld_win_dom.keys()
    }
    num_pairs_within = {
        annot: sum(num_win_dom[annot].values()) for annot in num_win_dom.keys()
    }

    ld_stats_between = {
        annot: np.sum(list(ld_btw_dom[annot].values()), axis=0)
        for annot in ld_btw_dom.keys()
    }
    num_pairs_between = {
        annot: sum(num_btw_dom[annot].values()) for annot in num_btw_dom.keys()
    }

    outputs = {
        "ld_within_domains": ld_stats_within,
        "num_pairs_within_domains": num_pairs_within,
        "ld_between_domains": ld_stats_between,
        "num_pairs_between_domains": num_pairs_between,
    }

    #    pickle.dump(outputs, open(f"parsed_data/{POP}.unphased.domains.bp", "wb+"))

    # bootstrap sigma_d^1 stats
    gene_sets = pickle.load(open("data/supp/bootstrap_gene_sets.500.bp", "rb"))
    # within domains
    bs_within = {}
    for annot in ld_stats_within.keys():
        bs_within[annot] = bootstrap_over_genes(ld_win_dom[annot], gene_sets)

    bs_between = {}
    for annot in ld_stats_between.keys():
        bs_between[annot] = bootstrap_over_genes(ld_btw_dom[annot], gene_sets)

    bs_outputs = {"bs_within": bs_within, "bs_between": bs_between}

    ## distances
    bin_edges = bin_edges = np.logspace(0, 7, 22)
    num_within = {
        annot: np.histogram(distances_within[annot], bin_edges) for annot in annots
    }
    num_between = {
        annot: np.histogram(distances_between[annot], bin_edges) for annot in annots
    }
    pickle.dump(
        {"bin_edges": bin_edges, "num_within": num_within, "num_between": num_between},
        open(f"parsed_data/{POP}.distances_within_between_domains.bp", "wb+"),
    )

#    pickle.dump(
#        bs_outputs,
#        open(
#            f"parsed_data/{POP}.unphased.domains.sigmad1.bootstrap_std_err.bp",
#            "wb+",
#        ),
#    )

# bootstrap mean D:
# within domains
# bs_within = {}
# for annot in ld_stats_within.keys():
#    bs_within[annot] = bootstrap_over_genes_meanD(
#        ld_win_dom[annot], num_win_dom[annot], gene_sets
#    )

# bs_between = {}
# for annot in ld_stats_between.keys():
#    bs_between[annot] = bootstrap_over_genes_meanD(
#        ld_btw_dom[annot], num_btw_dom[annot], gene_sets
#    )

# bs_outputs = {"bs_within": bs_within, "bs_between": bs_between}
