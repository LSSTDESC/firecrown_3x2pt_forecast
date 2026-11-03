#!/bin/bash
# Evaluate the likelihood once at the fiducial point in CosmoSIS, Cobaya and NumCosmo.
#
# Runs in the directories written by generate_configs.py and prints, for each
# framework, the log-likelihood (CosmoSIS, Cobaya) or -2 ln L (NumCosmo) it reports.
# Each run writes its output into <framework>/run_test.log.
#
# Usage: ./run_test.sh

cd "$(dirname "$0")" || exit 1

(cd cosmosis && cosmosis cosmosis_forecast_3x2pt.ini > run_test.log 2>&1)
echo "CosmoSIS (exit $?): $(grep -E '^Likelihood' cosmosis/run_test.log)"

(cd cobaya && cobaya-run -f cobaya_forecast_3x2pt.yaml > run_test.log 2>&1)
echo "Cobaya (exit $?): $(grep -E 'log-likelihood' cobaya/run_test.log)"

(cd numcosmo && numcosmo run test numcosmo_forecast_3x2pt.yaml > run_test.log 2>&1)
echo "NumCosmo (exit $?): $(grep -E 'm2lnL' numcosmo/run_test.log)"
