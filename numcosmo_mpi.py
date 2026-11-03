"""Run the numcosmo command line with the Firecrown types and models registered first.

Same command line as `numcosmo`: `python numcosmo_mpi.py run mcmc apes <experiment>.yaml ...`.

Under mpiexec, NumCosmo 0.28.0 starts the MPI worker loop while numcosmo_py.app is
imported (numcosmo_py/app/xcor/common.py calls Ncm.cfg_init() at module level). The
worker ranks then never reach numcosmo_py.app.loading, which registers the Firecrown
types and creates the models of <experiment>.builders.yaml, and they cannot deserialize
the work sent by rank 0. This script does both before importing numcosmo_py.app:

* it imports firecrown.connector.numcosmo, which registers FirecrownNumCosmoGaussCov
  and the other connector types, in every rank;
* in the worker ranks (MPI rank > 0) it creates the models of the builders file of the
  first .yaml argument. Rank 0 creates them in numcosmo_py.app.loading.
"""

import os
import sys
from pathlib import Path

import firecrown.connector.numcosmo  # noqa: F401
from numcosmo_py import Ncm

rank = int(os.environ.get("OMPI_COMM_WORLD_RANK", os.environ.get("PMI_RANK", "0")))
experiment = next((Path(a) for a in sys.argv[1:] if a.endswith(".yaml")), None)
if rank > 0 and experiment is not None:
    builders_file = experiment.with_suffix(".builders.yaml")
    if builders_file.exists():
        ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)  # as loading.py in 0.28.0
        model_builders = ser.dict_str_from_yaml_file(builders_file.absolute().as_posix())
        for name in model_builders.keys():  # noqa: SIM118
            model_builders.get(name).create()

from numcosmo_py.app import app  # noqa: E402

app()
