#!/bin/bash
# Evaluate the likelihood once at the fiducial point in CosmoSIS, Cobaya and NumCosmo.
#
# Runs in the directories written by generate_configs.py and prints, for each
# framework, the log-likelihood (CosmoSIS, Cobaya) or -2 ln L (NumCosmo) it reports.
# Each run writes its output into <framework>/run_test.log. Cobaya's evaluate sampler
# would draw its point from the ref distribution of each free parameter, so the test
# writes cobaya/cobaya_forecast_3x2pt_test.yaml, which sets the point to the ref means
# (the fiducial values) with the evaluate sampler's override.
#
# Usage: ./run_test.sh

cd "$(dirname "$0")" || exit 1

(cd cosmosis && cosmosis cosmosis_forecast_3x2pt.ini > run_test.log 2>&1)
echo "CosmoSIS (exit $?): $(grep -E '^Likelihood' cosmosis/run_test.log)"

(cd cobaya && python - <<'EOF' && cobaya-run -f cobaya_forecast_3x2pt_test.yaml > run_test.log 2>&1)
import re
text = open("cobaya_forecast_3x2pt.yaml").read()
assert "sampler:\n  evaluate: null\n" in text
means = re.findall(r"\n  (\w+):\n    ref:\n      dist: norm\n      loc: (\S+)\n", text)
assert means
override = "".join(f"      {name}: {loc}\n" for name, loc in means)
text = text.replace("sampler:\n  evaluate: null\n", "sampler:\n  evaluate:\n    override:\n" + override)
open("cobaya_forecast_3x2pt_test.yaml", "w").write(text)
EOF
echo "Cobaya (exit $?): $(grep -E 'log-likelihood' cobaya/run_test.log)"

(cd numcosmo && numcosmo run test numcosmo_forecast_3x2pt.yaml > run_test.log 2>&1)
echo "NumCosmo (exit $?): $(grep -E 'm2lnL' numcosmo/run_test.log)"
