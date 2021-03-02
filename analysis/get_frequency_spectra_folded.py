# Get the frequency spectra for each mutation type, for all mutations
# within a gene, for mutations falling in annotated domains, and for
# mutations falling outside of annotated domains.

# Ignore ancestral information, so keeping all SNPs, and return folded SFS

import numpy as np
import gzip
import moments
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


def get_annotations(annot_file):
    annots = {}  # pos: alt: (gene, csq)
    for line in open(annot_file, "r"):
        chrom, pos, ref, alt, gene, csq = line.split()
        annots.setdefault(pos, {})
        annots[pos][alt] = (gene, csq)
    return annots


def to_haplotypes(gts):
    h = np.zeros(2 * len(gts))
    for ii, g in enumerate(gts):
        h[2 * ii] = int(g[0])
        h[2 * ii + 1] = int(g[2])
    return h


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
    kept = 0
    discarded = 0
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
            gts = l.split()[9:]
            if np.all([g == "1|1" for g in gts]) or np.all([g == "0|0" for g in gts]):
                # fixed for ref or alt
                continue
            gene, csq = annot_dict[pos][alt]
            positions.append(int(pos))
            annotations.append(csq)
            genes.append(gene)
            G.append(to_haplotypes(gts))

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


def compute_frequency_spectra(pos, G, domains):
    domain_ls = []
    domain_rs = []
    for gene, intervals in domains[chrom].items():
        for interval in intervals:
            assert interval[1] >= interval[0]
            domain_ls.append(interval[0])
            domain_rs.append(interval[1])
    domain_ls = np.array(domain_ls)
    domain_rs = np.array(domain_rs)
    in_domain = [np.any(np.logical_and(domain_ls <= p, p <= domain_rs)) for p in pos]
    G_in = G.compress(in_domain, axis=0)
    G_out = G.compress(np.logical_not(in_domain), axis=0)
    fs_all = np.zeros(len(G[0]) + 1)
    fs_in = np.zeros(len(G[0]) + 1)
    fs_out = np.zeros(len(G[0]) + 1)
    freqs, counts = np.unique(G.sum(axis=1), return_counts=True)
    fs_all[freqs.astype(int)] = counts
    freqs, counts = np.unique(G_in.sum(axis=1), return_counts=True)
    fs_in[freqs.astype(int)] = counts
    freqs, counts = np.unique(G_out.sum(axis=1), return_counts=True)
    fs_out[freqs.astype(int)] = counts
    return fs_all, fs_in, fs_out


if __name__ == "__main__":
    POP = sys.argv[1]

    domains = pickle.load(open("data/supp/domain_dict.bp", "rb"))

    annots = ["synonymous", "missense", "loss_of_function"]
    spectra = {annot: {"all": [], "in": [], "out": []} for annot in annots}

    for chrom in range(1, 23):
        print("processing chromosome", chrom)
        positions, annotations, genes, G = load_genotype_matrix(chrom, POP)
        for annot in annots:
            pos_sub, genes_sub, G_sub = subset_annotation(
                positions, annotations, genes, G, annot=annot
            )
            chrom_spectra = compute_frequency_spectra(pos_sub, G_sub, domains)
            spectra[annot]["all"].append(chrom_spectra[0])
            spectra[annot]["in"].append(chrom_spectra[1])
            spectra[annot]["out"].append(chrom_spectra[2])

    compiled_spectra = {}
    for annot in annots:
        compiled_spectra[annot] = {}
        for category in ["all", "in", "out"]:
            fs = moments.Spectrum(np.sum(spectra[annot][category], axis=0))
            fs.fold()
            compiled_spectra[annot][category] = fs

    for annot in annots:
        assert np.all(compiled_spectra[annot]["all"] == compiled_spectra[annot]["in"] + compiled_spectra[annot]["out"])

    pickle.dump(compiled_spectra, open(f"parsed_data/{POP}.frequency_spectra.folded.bp", "wb+"))
