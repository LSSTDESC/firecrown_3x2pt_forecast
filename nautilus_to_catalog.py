#!/usr/bin/env python
"""Write the nautilus samples as a weighted NumCosmo catalog.

Reads cosmosis/output/nautilus.txt and writes cosmosis/output/nautilus.mcmc.fits, an
NcmMSetCatalog on the model set of the NumCosmo experiment (whose free parameters are the
8 CosmoSIS free parameters) with one chain and one row per nautilus point. Each row holds
-2 ln L = -2 (post - prior), the weight exp(log_weight - max log_weight) and the
parameters, so `numcosmo catalog analyze` and `numcosmo catalog plot-corner` read it as
they read the APES catalog. NumCosmo computes no autocorrelation time for a weighted
catalog; the nautilus points are independent weighted draws. NumCosmo 0.30 aborts on a
weighted catalog wherever it estimates the autocorrelation time: `catalog analyze` and
the automatic thinning of `catalog plot-corner`, so plot-corner needs `--thin 1`. Its
`--burnin` applies to every catalog; the first 18 nautilus rows have zero weight.

Usage: python nautilus_to_catalog.py
"""

import os
from pathlib import Path

import numpy as np
from numcosmo_py import Ncm
from numcosmo_py.app.loading import LoadExperiment

HERE = Path(__file__).resolve().parent
SAMPLES = HERE / "cosmosis" / "output" / "nautilus.txt"
OUTPUT = HERE / "cosmosis" / "output" / "nautilus.mcmc.fits"

with SAMPLES.open(encoding="utf-8") as f:
    header = f.readline().lstrip("#").split()
data = np.loadtxt(SAMPLES)

# The experiment loads its likelihood factory relative to the working directory.
os.chdir(HERE / "numcosmo")
mset = LoadExperiment(experiment=HERE / "numcosmo" / "numcosmo_forecast_3x2pt.yaml").mset
mset.prepare_fparam_map()

# CosmoSIS columns are <section>--<name> with lower-case names.
cosmosis_column = {name.split("--")[-1].lower(): i for i, name in enumerate(header)}
columns = [cosmosis_column[mset.fparam_name(i).lower()] for i in range(mset.fparams_len())]
m2lnL = -2.0 * (data[:, header.index("post")] - data[:, header.index("prior")])
log_weight = data[:, header.index("log_weight")]
weight = np.exp(log_weight - log_weight.max())

OUTPUT.unlink(missing_ok=True)
storage = Ncm.MSetCatalogStorageFits.new(OUTPUT.as_posix(), False)
desc = Ncm.MSetCatalogDesc.new(
    mset, 1, True, [Ncm.MSET_CATALOG_M2LNL_COLNAME], [Ncm.MSET_CATALOG_M2LNL_SYMBOL], 0, "nautilus"
)
storage.set_desc(desc)
mcat = Ncm.MSetCatalog.new_from_storage(storage)
for row, m2lnL_i, weight_i in zip(data[:, columns], m2lnL, weight):
    mcat.add_from_vector_array(Ncm.Vector.new_array(row.tolist()), [float(m2lnL_i), float(weight_i)])
mcat.sync(True)
print(f"{len(data)} nautilus points written to {OUTPUT}")
