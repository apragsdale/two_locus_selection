# for each protein (ensembl ID), get list of regions based on prot2hg file
from collections import defaultdict

data = {}
for i in range(1, 23):
    data[i] = defaultdict(list)

c = 0
Ls = []

chroms_to_keep = [str(i) for i in range(1, 23)]

header = 0
for line in open("prot2hg_1938_112019.csv", "r"):
    if header == 0:
        chr_col = line.split(";").index('"chrom"')
        start_col = line.split(";").index('"chr_start"')
        end_col = line.split(";").index('"chr_end"')
        ensembl_col = line.split(";").index('"ensembl"')
        type_col = line.split(";").index('"type"')
        header = 1
    else:
        l = line.split(";")
        t = l[type_col]
        if t != '"Region"':
            continue
        chrom = l[chr_col]
        chrom = chrom.split('"')[1][3:]
        if chrom not in chroms_to_keep:
            continue
        start = int(l[start_col])
        end = int(l[end_col])
        ID = l[ensembl_col]
        ID = ID.split('"')[1]
        IDs = ID.split(",")
        for ID in IDs:
            data[int(chrom)][ID].append((start, end))
        c += 1

print("found", c, "domains")

import pickle
pickle.dump(data, open("domain_dict.bp", "wb+"))
