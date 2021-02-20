import gzip

coding_annotations = [
    "frameshift_variant",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "start_lost",
    "stop_gained",
    "stop_lost",
    "transcript_ablation",
    "missense",
    "synonymous",
]

lof_annotations = coding_annotations[:-2]

SNPs = ["A", "C", "G", "T"]

for chrom in range(1, 23):
    with gzip.open(
        f"../../../../../Data/1000G/annotations/ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5.20130502.sites.annotation.vcf.gz",
        "rb"
    ) as fin:
        with open(f"variants/coding_variant_consequences.chr{chrom}.txt", "w+") as fout:
            for line in fin:
                l = line.decode()
                if not l.startswith("#"):
                    coding = False
                    for annot in coding_annotations:
                        if annot in l:
                            coding = True
                    if coding:
                        c, p, _, ref, alts, _, _, info = l.split()
                        if ref not in SNPs:
                            continue
                        csqs = info.split("CSQ=")[1].split(";")[0].split(",")
                        alts = alts.split(",")
                        for alt in alts:
                            if alt not in SNPs:
                                continue
                            found = False
                            for csq in csqs:
                                allele, gene, _, _, annot = csq.split("|")[:5]
                                if allele == alt:
                                    found = True
                                    fout.write(f"{c}\t{p}\t{ref}\t{alt}\t{gene}\t{annot}\n")
                            if not found:
                                raise ValueError("uh oh")
    print("finished with chromosome", chrom)
