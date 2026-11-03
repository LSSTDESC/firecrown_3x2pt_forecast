#!/usr/bin/env python
"""Write the posterior covariance of the APES chain as a Cobaya covmat file.

Reads numcosmo/apes.mcmc.fits, drops the first BURNIN iterations (BURNIN x NWALKERS
rows; the catalog stores the walkers of one iteration in consecutive rows), and writes
the covariance of the free parameters to cobaya/apes.covmat in Cobaya's format: a header
line `# name1 name2 ...` with the Cobaya parameter names, then the matrix.

Usage: python apes_covmat.py [BURNIN]   (default 12: where the variance of -2 ln L
                                         across the ensemble stops decreasing)
"""

import sys
from pathlib import Path

import numpy as np
from astropy.io import fits

HERE = Path(__file__).resolve().parent
CATALOG = HERE / "numcosmo" / "apes.mcmc.fits"
OUTPUT = HERE / "cobaya" / "apes.covmat"
NWALKERS = 900
BURNIN = int(sys.argv[1]) if len(sys.argv) > 1 else 12

with fits.open(CATALOG) as hdus:
    table = hdus[1].data
    columns = [c for c in hdus[1].columns.names if c != "NcmFit:m2lnL"]
    samples = np.column_stack([table[c] for c in columns])[BURNIN * NWALKERS :]

names = [c.split(":", 1)[1] for c in columns]
covariance = np.cov(samples, rowvar=False)
np.savetxt(OUTPUT, covariance, header=" ".join(names), comments="# ")
print(f"{samples.shape[0]} samples after {BURNIN} iterations of burn-in -> {OUTPUT}")
for name, sigma in zip(names, np.sqrt(np.diag(covariance))):
    print(f"  {name:12s} sigma = {sigma:.6g}")
