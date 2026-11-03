#!/bin/bash
# Run the NumCosmo APES sampler on the 3x2pt forecast with MPI.
#
# Runs numcosmo run mcmc apes on numcosmo/numcosmo_forecast_3x2pt.yaml (written by
# generate_configs.py) under mpiexec with NPROC processes: rank 0 runs APES and the
# other NPROC - 1 ranks evaluate the likelihood, one OpenMP thread each. It needs the
# MPI build of NumCosmo (conda-forge build string mpi_openmpi or mpi_mpich). The
# command line is started through numcosmo_mpi.py, which registers the Firecrown types
# in every rank before NumCosmo starts the MPI workers. The catalog is written to
# numcosmo/apes.mcmc.fits and the log to numcosmo/apes.log.
#
# When numcosmo/apes.mcmc.fits exists, APES continues it up to NSAMPLES iterations, with
# the random number generator state saved in the catalog; SEED is then not passed, since
# NumCosmo refuses to replace that state. Delete the catalog to start a new chain.
#
# Settings, each overridable from the environment:
#   NPROC     MPI processes (default 11: 1 master and 10 likelihood workers)
#   NWALKERS  APES walkers (default 0: NumCosmo chooses the number)
#   NSAMPLES  samples per walker, i.e. iterations of the ensemble (default 50); the
#             catalog has NWALKERS x NSAMPLES points
#   SEED      seed of the sampler and of the initial points of a new chain (default 1)
#
# Usage: ./run_apes.sh

cd "$(dirname "$0")/numcosmo" || exit 1

NPROC=${NPROC:-11}
NWALKERS=${NWALKERS:-0}
NSAMPLES=${NSAMPLES:-50}
SEED=${SEED:-1}

export OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=1 OPENBLAS_NUM_THREADS=1 MPLBACKEND=Agg

SEED_ARGS=(--seed "$SEED")
if [ -e apes.mcmc.fits ]; then
    SEED_ARGS=()
fi

# Open MPI counts physical cores as slots; --use-hwthread-cpus counts hardware threads,
# so NPROC may exceed the number of cores. --bind-to none lets the operating system
# place the processes. --mca btl self,sm restricts communication to shared memory, the
# transport between processes on one machine.
mpiexec -n "$NPROC" --use-hwthread-cpus --bind-to none --mca btl self,sm \
    python ../numcosmo_mpi.py run mcmc apes numcosmo_forecast_3x2pt.yaml \
    --parallel mpi \
    --nwalkers "$NWALKERS" \
    --nsamples "$NSAMPLES" \
    "${SEED_ARGS[@]}" \
    --output apes.fits \
    --log-file apes.log
