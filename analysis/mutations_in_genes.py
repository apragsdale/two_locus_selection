import numpy as np
import gzip
import moments.LD
import sys
import pickle
from collections import defaultdict
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


lof_annotations = [
    "frameshift_variant",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "start_lost",
    "stop_gained",
    "stop_lost",
    "transcript_ablation",
]


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


if __name__ == "__main__":
    chroms = range(1, 23)
    annots = ["synonymous", "missense", "loss_of_function"]
    pops = {
        "AFR": ["YRI", "ESN", "GWD", "LWK", "MSL"],
        "EAS": ["CHB", "CHS", "CDX", "KHV", "JPT"],
        "EUR": ["CEU", "GBR", "TSI", "IBS", "FIN"],
    }

    # stores set for each gene with positions of that mutation type
    mutations_by_gene = {annot: defaultdict(set) for annot in annots}
    for _, v in pops.items():
        for POP in v:
            print(POP)
            for chrom in chroms:
                print(" ", chrom)
                positions, annotations, genes, G = load_genotype_matrix(chrom, POP)
                for annot in annots:
                    pos_sub, genes_sub, G_sub = subset_annotation(
                        positions, annotations, genes, G, annot=annot
                    )
                    gene_set = set(genes_sub)
                    for gene in gene_set:
                        to_keep = genes_sub == gene
                        pos_gene = pos_sub.compress(to_keep)
                        for pos in pos_gene:
                            mutations_by_gene[annot][gene].add(pos)

    all_genes = set.union(
        set(mutations_by_gene["synonymous"].keys()),
        set(mutations_by_gene["missense"].keys()),
    )

    mutation_counts = {}
    # [num_syn, num_mis]
    for gene in all_genes:
        mutation_counts[gene] = [
            len(mutations_by_gene["synonymous"][gene]),
            len(mutations_by_gene["missense"][gene]),
        ]

    count_pairs = np.array(list(mutation_counts.values()))

    corr = np.corrcoef(count_pairs.T)[0, 1]
    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(
        count_pairs[:, 0], count_pairs[:, 1]
    )

    jitter = np.random.rand(np.prod(count_pairs.shape)).reshape(count_pairs.shape) - 0.5

    fig = plt.figure(54329, figsize=(5, 4))
    fig.clf()
    ax = plt.subplot(111)
    ax.scatter(
        count_pairs[:, 0] + jitter[:, 0],
        count_pairs[:, 1] + jitter[:, 1],
        marker=".",
        s=1,
        alpha=1,
    )
    xlim = ax.get_xlim()
    xx = np.linspace(xlim[0], xlim[1])
    yy = intercept + slope * xx
    ax.plot(xx, yy, "k--", label=f"$y={slope:.2f}x +{intercept:.2f}$\n$r={corr:.2f}$")
    ax.set_xlim(-1, 121)
    ax.set_ylim(-1, 241)
    # ax.set_xscale("log")
    # ax.set_yscale("log")
    ax.set_xlabel("Number of synonymous mutations")
    ax.set_ylabel("Number of missense mutations")
    ax.legend()
    ax.set_title("Numbers of mutations within the same gene")
    fig.tight_layout()
    plt.savefig("../figures/mutations_in_genes.pdf")
