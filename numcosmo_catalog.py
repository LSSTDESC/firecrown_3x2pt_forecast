#!/usr/bin/env python
"""Run a `numcosmo catalog` command on catalogs of the 3x2pt forecast.

Creates the models of numcosmo/numcosmo_forecast_3x2pt.builders.yaml, then runs
`numcosmo catalog` with the arguments given. A catalog stores the model set, but the
model types FirecrownCclCosmo and FirecrownTwoPoint exist only after their builders are
created, and NumCosmo 0.28.0 does not create them when it reads a catalog.

Usage:
    python numcosmo_catalog.py analyze numcosmo/apes.mcmc.fits --burnin 10
    python numcosmo_catalog.py plot-corner numcosmo/apes.mcmc.fits --burnin 10 \
        --no-show --output numcosmo/apes
The second command writes numcosmo/apes.corner.pdf: plot-corner replaces the suffix of
--output by .corner.pdf.
"""

import sys
from pathlib import Path

from numcosmo_py import Ncm

HERE = Path(__file__).resolve().parent
BUILDERS_FILE = HERE / "numcosmo" / "numcosmo_forecast_3x2pt.builders.yaml"

ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
model_builders = ser.dict_str_from_yaml_file(BUILDERS_FILE.as_posix())
for name in model_builders.keys():  # noqa: SIM118
    model_builders.get(name).create()

from numcosmo_py.app import app  # noqa: E402

sys.argv = ["numcosmo", "catalog", *sys.argv[1:]]
app()
