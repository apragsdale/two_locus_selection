import pickle
import moments
import numpy as np

pops = ["ESN", "GWD", "LWK", "MSL", "YRI", "CEU", "FIN", "GBR", "IBS", "TSI", "CDX", "CHB", "CHS", "JPT", "KHV"]

with open("tajimas_d.txt", "w+") as fout:
    fout.write("Population;Mutation type;Region;Tajima's D\n")
    for pop in pops:
        spectra = pickle.load(open(f"../analysis/parsed_data/{pop}.frequency_spectra.bp", "rb"))
        Ds = {}
        for annot in ["synonymous", "missense", "loss_of_function"]:
            Ds[annot] = {}
            for c in ["all", "in", "out"]:
                Ds[annot][c] = moments.Spectrum(spectra[annot][c]).Tajima_D()
        fout.write(f"{pop};Synonymous;All;{Ds['synonymous']['all']:.3f}\n")
        fout.write(f";;In domain;{Ds['synonymous']['in']:.3f}\n")
        fout.write(f";;Not in domain;{Ds['synonymous']['out']:.3f}\n")
        fout.write(f";Missense;All;{Ds['missense']['all']:.3f}\n")
        fout.write(f";;In domain;{Ds['missense']['in']:.3f}\n")
        fout.write(f";;Not in domain;{Ds['missense']['out']:.3f}\n")
        fout.write(f";Loss of function;All;{Ds['loss_of_function']['all']:.3f}\n")
        fout.write(f";;In domain;{Ds['loss_of_function']['in']:.3f}\n")
        fout.write(f";;Not in domain;{Ds['loss_of_function']['out']:.3f}\n")

