#!/bin/bash
# Run the CosmoSIS nautilus sampler on the 3x2pt forecast with MPI.
#
# Runs cosmosis --mpi on cosmosis/cosmosis_forecast_3x2pt.ini (written by
# generate_configs.py) under mpiexec with NPROC processes, one OpenMP thread each, with
# the sampler set to nautilus on the command line. The samples are written to
# cosmosis/output/nautilus.txt and the log to cosmosis/nautilus.log.
#
# Settings, each overridable from the environment:
#   NPROC   MPI processes (default 11)
#   NLIVE   nautilus live points (default 1500, the value of the minimal_mcp_pipelines
#           Perlmutter job)
#
# Usage: ./run_nautilus.sh

cd "$(dirname "$0")/cosmosis" || exit 1

NPROC=${NPROC:-11}
NLIVE=${NLIVE:-1500}

export OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=1 OPENBLAS_NUM_THREADS=1 MPLBACKEND=Agg

# Same mpiexec options as run_apes.sh: hardware threads as slots, no binding, shared
# memory transport.
mpiexec -n "$NPROC" --use-hwthread-cpus --bind-to none --mca btl self,sm \
    cosmosis --mpi cosmosis_forecast_3x2pt.ini \
    -p runtime.sampler=nautilus \
    nautilus.n_live="$NLIVE" \
    nautilus.verbose=T \
    output.filename=output/nautilus.txt \
    > nautilus.log 2>&1
