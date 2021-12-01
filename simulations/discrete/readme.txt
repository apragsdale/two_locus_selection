The scripts in this directory compare expectations from moments.TwoLocus to
discrete simulations. By "discrete", we mean two loci separated by a given
recombination distance and no inference from linked selected alleles (see the
simulations in ../fwdpy11/ for those simulation comparisons).

This scenario is primarily for verification: single-population of constant
size, given selection and dominance or epistasis coefficients. For each
scenario, we focus on residuals between expectations and simulations of the
full two-locus haplotype frequency spectrum.

Scenarios:

gamma = -0.1, -2, -20

(with N = 5000, s = gamma / 2 / N, so s = -1e-5, -2e-4, -2e-3)

rho = 0.1, 1, 20

(with N = 5000, r = rho / 4 / N, so r = 5e-6, 5e-5, 1e-3)

This gives 9 combinations of gamma x rho (assuming selection is
equal at the left and right loci).

And then consider either epistasis or dominance:

e = -0.5, 0, 0.5 (with h = 0.5)

(or model is sAB = (sA + sB) * (1 + e), so e > 0 is synergistic
and e < 0 is diminishing returns)

h = 0.1 (with e = 0)

This results in 45 different scenarios.
