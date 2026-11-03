#!/bin/bash
# Run the Cobaya mcmc sampler on the 3x2pt forecast with MPI.
#
# Writes cobaya/cobaya_forecast_3x2pt_mcmc.yaml from cobaya/cobaya_forecast_3x2pt.yaml
# (written by generate_configs.py) with the sampler replaced by mcmc, then runs it under
# mpiexec with NPROC processes: Cobaya runs one chain per process and stops when the
# Gelman-Rubin R - 1 of the means across chains is below Rminus1_stop (0.01, Cobaya's
# default). The free parameters start from the fiducial point (written by
# generate_configs.py). When cobaya/apes.covmat exists (written by apes_covmat.py from
# the APES chain), it is the initial proposal covariance; otherwise the proposal widths
# written by generate_configs.py are used. The chains are written to
# cobaya/mcmc/chain.*.txt and the log to cobaya/mcmc.log (appended to).
#
# When cobaya/mcmc/ exists, Cobaya continues those chains from their checkpoint
# (cobaya-run --resume); delete cobaya/mcmc/ to start new chains. A run that was killed
# leaves cobaya/mcmc/chain.input.yaml.locked behind; delete it before resuming when no
# Cobaya process is running.
#
# Settings, each overridable from the environment:
#   NPROC  MPI processes, i.e. chains (default 10)
#
# Usage: ./run_cobaya.sh

cd "$(dirname "$0")/cobaya" || exit 1

NPROC=${NPROC:-10}

python - <<'EOF' || exit 1
from pathlib import Path
text = open("cobaya_forecast_3x2pt.yaml").read()
assert "sampler:\n  evaluate: null\n" in text and "\noutput: output\n" in text
mcmc = "  mcmc:\n    covmat: apes.covmat\n" if Path("apes.covmat").exists() else "  mcmc: {}\n"
text = text.replace("sampler:\n  evaluate: null\n", "sampler:\n" + mcmc)
text = text.replace("\noutput: output\n", "\noutput: mcmc/chain\n")
open("cobaya_forecast_3x2pt_mcmc.yaml", "w").write(text)
EOF

export OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=1 OPENBLAS_NUM_THREADS=1 MPLBACKEND=Agg

MODE=--force
if [ -d mcmc ]; then
    MODE=--resume
fi

# Same mpiexec options as run_apes.sh: hardware threads as slots, no binding, shared
# memory transport.
mpiexec -n "$NPROC" --use-hwthread-cpus --bind-to none --mca btl self,sm \
    cobaya-run "$MODE" cobaya_forecast_3x2pt_mcmc.yaml >> mcmc.log 2>&1
