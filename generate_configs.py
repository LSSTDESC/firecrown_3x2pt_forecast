#!/usr/bin/env python
"""Generate CosmoSIS, Cobaya and NumCosmo configurations for the 3x2pt forecast.

For each framework, deletes the directory <framework>/ and writes into a new one a copy
of the SACC file, a copy of the likelihood factory forecast_3x2pt.py and the framework
configuration. The
configurations evaluate the likelihood once at the fiducial point (CosmoSIS sampler
"test", Cobaya sampler "evaluate", NumCosmo "numcosmo run test").

The SACC file is forecast_3x2pt_srd_y1_noisy.sacc, or forecast_3x2pt_srd_y1.sacc with
--noiseless; generate_data.py writes both. The cosmology is computed by CCL inside the
likelihood (FrameworkCosmology.NONE), so the frameworks sample the CCL parameters of
the model firecrown_ccl_cosmo and run no Boltzmann code.

Free parameters, with uniform priors: Omega_c, sigma8, ia_bias and lens{i}_bias. Every
other parameter is fixed at its fiducial value in generate_data.FIDUCIAL.

Cobaya starts from the information APES starts from: for each free parameter, the
starting point is drawn from a Gaussian centered on the fiducial value with standard
deviation equal to the parameter scale (1% of the fiducial value), and the scale is the
initial proposal width.

Three changes to what Firecrown 1.16.1 writes:

* its generator writes values with 3 significant digits
  (firecrown.app.analysis._cosmosis.format_float). This script replaces format_float by
  a function that writes the shortest exact decimal representation, so that the fixed
  parameters keep their fiducial values;
* its .ini points to the connector inside the Firecrown installation. This script
  writes cosmosis/firecrown_cosmosis.py, which imports the connector of the Firecrown
  installed in the running Python environment, and points the .ini to it, as the
  generator on Firecrown master does;
* its NumCosmo experiment file sets <name>-fit: false for every parameter of the model
  set, which overrides the free parameters declared in the builders file. This script
  sets <name>-fit: true for the free parameters.

Usage: python generate_configs.py [--noiseless] [--only cosmosis|cobaya|numcosmo]
--only regenerates one framework directory and leaves the others untouched.
"""

import shutil
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import ClassVar

import firecrown.app.analysis._cosmosis as _cosmosis
from firecrown.app.analysis import (
    AnalysisBuilder,
    FrameworkCosmology,
    Frameworks,
    Model,
    Parameter,
)
from firecrown.app.sacc import SaccFormat
from firecrown.likelihood import NamedParameters

from generate_data import FIDUCIAL, OUTPUT_NOISY_SACC, OUTPUT_SACC

_cosmosis.format_float = lambda value: repr(float(value))

HERE = Path(__file__).resolve().parent
SACC_FILE = OUTPUT_SACC if "--noiseless" in sys.argv[1:] else OUTPUT_NOISY_SACC
FACTORY_FILE = HERE / "forecast_3x2pt.py"

# (name, LaTeX symbol, lower bound, upper bound, free); the value is FIDUCIAL[name].
CCL_PARAMETERS = [
    ("Omega_c", r"\Omega_c", 0.06, 0.46, True),
    ("sigma8", r"\sigma_8", 0.3, 1.2, True),
    ("Omega_b", r"\Omega_b", 0.03, 0.07, False),
    ("h", "h", 0.5, 0.8, False),
    ("n_s", "n_s", 0.85, 1.10, False),
    ("Omega_k", r"\Omega_k", -0.2, 0.2, False),
    ("w0", "w_0", -3.0, 0.0, False),
    ("wa", "w_a", -1.0, 1.0, False),
    ("Neff", r"N_{\rm eff}", 2.0, 5.0, False),
    ("m_nu", r"m_\nu", 0.0, 5.0, False),
    ("T_CMB", r"T_{\rm CMB}", 2.7, 2.8, False),
]
FIRECROWN_PARAMETERS = [
    ("ia_bias", r"A_{\rm IA}", -5.0, 5.0, True),
    ("alphaz", r"\alpha_z", -5.0, 5.0, False),
    ("z_piv", r"z_{\rm piv}", 0.0, 1.0, False),
]
for i in range(5):
    FIRECROWN_PARAMETERS += [
        (f"lens{i}_bias", f"b_{i}", 0.8, 2.5, True),
        (f"lens{i}_mag_bias", f"s_{i}", -0.01, 0.01, False),
        (f"lens{i}_delta_z", rf"\delta z^{{\rm lens}}_{i}", -0.05, 0.05, False),
        (f"lens{i}_alphag", rf"\alpha_{{g,{i}}}", -5.0, 5.0, False),
        (f"lens{i}_alphaz", rf"\alpha_{{z,{i}}}", -5.0, 5.0, False),
        (f"lens{i}_z_piv", rf"z_{{\rm piv,{i}}}", 0.0, 1.0, False),
        (f"src{i}_delta_z", rf"\delta z^{{\rm src}}_{i}", -0.05, 0.05, False),
    ]


@dataclass
class Forecast3x2pt(AnalysisBuilder):
    """Analysis builder: copies SACC_FILE and FACTORY_FILE and declares the parameters."""

    description: ClassVar[str] = "3x2pt forecast, LSST DESC SRD Y1 noise"

    def generate_sacc(self, output_path: Path) -> Path:
        """Copy SACC_FILE into output_path."""
        return Path(shutil.copy(SACC_FILE, output_path / SACC_FILE.name))

    def generate_factory(self, output_path: Path, sacc: Path) -> Path:
        """Copy FACTORY_FILE into output_path."""
        return Path(shutil.copy(FACTORY_FILE, output_path / FACTORY_FILE.name))

    def get_build_parameters(self, sacc_path: Path) -> NamedParameters:
        """Pass the SACC file to build_likelihood."""
        return NamedParameters({"sacc_file": self.get_sacc_file(sacc_path)})

    def get_models(self) -> list[Model]:
        """Return the CCL parameters and the Firecrown systematics parameters."""
        return [
            Model(
                name="firecrown_ccl_cosmo",
                description="CCL cosmological parameters",
                parameters=[
                    Parameter.from_tuple(n, s, lo, hi, FIDUCIAL[n], free)
                    for n, s, lo, hi, free in CCL_PARAMETERS
                ],
            ),
            Model(
                name="firecrown_two_point",
                description="Galaxy bias, intrinsic alignment, magnification and photo-z",
                parameters=[
                    Parameter.from_tuple(n, s, lo, hi, FIDUCIAL[n], free)
                    for n, s, lo, hi, free in FIRECROWN_PARAMETERS
                ],
            ),
        ]

    def required_cosmology(self) -> FrameworkCosmology:
        """The likelihood computes the cosmology with CCL."""
        return FrameworkCosmology.NONE


if __name__ == "__main__":
    frameworks = list(Frameworks)
    if "--only" in sys.argv:
        frameworks = [Frameworks(sys.argv[sys.argv.index("--only") + 1])]
    free_parameters = [
        Parameter.from_tuple(n, s, lo, hi, FIDUCIAL[n], free)
        for n, s, lo, hi, free in CCL_PARAMETERS + FIRECROWN_PARAMETERS
        if free
    ]

    for framework in frameworks:
        shutil.rmtree(HERE / framework.value, ignore_errors=True)
        Forecast3x2pt(
            output_path=HERE / framework.value,
            prefix="forecast_3x2pt",
            target_framework=framework,
            sacc_format=SaccFormat.FITS,
            use_absolute_path=False,
            quiet=True,
        )

    if Frameworks.COSMOSIS in frameworks:
        # Point the CosmoSIS .ini to the connector shim instead of the installation path.
        cosmosis_dir = HERE / Frameworks.COSMOSIS.value
        (cosmosis_dir / "firecrown_cosmosis.py").write_text(
            "from firecrown.connector.cosmosis.likelihood import cleanup, execute, setup\n",
            encoding="utf-8",
        )
        ini_path = cosmosis_dir / "cosmosis_forecast_3x2pt.ini"
        lines = ini_path.read_text(encoding="utf-8").splitlines()
        lines = [
            "file = firecrown_cosmosis.py"
            if line.startswith("file =") and line.endswith("connector/cosmosis/likelihood.py")
            else line
            for line in lines
        ]
        ini_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    if Frameworks.COBAYA in frameworks:
        # Start Cobaya from the information APES starts from: the starting point is drawn
        # from a Gaussian centered on the fiducial value with standard deviation scale,
        # and scale is also the initial proposal width.
        yaml_path = HERE / Frameworks.COBAYA.value / "cobaya_forecast_3x2pt.yaml"
        text = yaml_path.read_text(encoding="utf-8")
        for p in free_parameters:
            old = f"  {p.name}:\n    ref: {p.default_value!r}\n"
            assert text.count(old) == 1, p.name
            new = (
                f"  {p.name}:\n    ref:\n      dist: norm\n      loc: {p.default_value!r}\n"
                f"      scale: {p.scale!r}\n    proposal: {p.scale!r}\n"
            )
            text = text.replace(old, new)
        yaml_path.write_text(text, encoding="utf-8")

    if Frameworks.NUMCOSMO in frameworks:
        # Mark the free parameters free in the NumCosmo model set. NumCosmo property names
        # replace "_" by "-".
        experiment_path = HERE / Frameworks.NUMCOSMO.value / "numcosmo_forecast_3x2pt.yaml"
        text = experiment_path.read_text(encoding="utf-8")
        for p in free_parameters:
            flag = f"{p.name.replace('_', '-')}-fit:"
            assert text.count(f" {flag} false\n") == 1, flag
            text = text.replace(f" {flag} false\n", f" {flag} true\n")
        experiment_path.write_text(text, encoding="utf-8")
