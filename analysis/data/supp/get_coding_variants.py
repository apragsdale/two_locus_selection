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

for chrom in range(1, 23):
    with gzip.open(
        f"../../../../../Data/1000G/annotations/ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5.20130502.sites.annotation.vcf.gz",
        "rb"
    ) as fin:
        with open(f"variants/coding_variants.chr{chrom}.txt", "w+") as fout:
            for line in fin:
                l = line.decode()
                if not l.startswith("#"):
                    coding = False
                    for annot in coding_annotations:
                        if annot in l:
                            coding = True
                    if coding:
                        c, p = l.split()[:2]
                        fout.write(f"{c}\t{p}\n")
    print("finished with chromosome", chrom)
