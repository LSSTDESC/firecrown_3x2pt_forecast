#!/usr/bin/env python
"""Generate the data vectors and covariance of the 3x2pt forecast.

Writes two SACC files with the same tracers, spectra, bandpower windows and covariance:

* forecast_3x2pt_srd_y1.sacc: the mean is the theory vector of forecast_3x2pt.py at
  the fiducial point;
* forecast_3x2pt_srd_y1_noisy.sacc: the mean plus L z, where C = L L^T is the Cholesky
  factorization of the covariance and z is standard normal, drawn after
  np.random.seed(NOISE_SEED).

The tracers, their n(z), the 27 spectra and the 20 top-hat bandpower windows
(20 <= ell <= 2000) are read from INPUT_SACC. Its data vector and covariance are not
used.

The covariance is the Gaussian covariance of the bandpowers on a fraction F_SKY of
the sky,

    Cov[C^ab_L, C^cd_L'] = sum_ell w_L(ell) w_L'(ell)
                           [Ct^ac Ct^bd + Ct^ad Ct^bc] / ((2 ell + 1) F_SKY),

where w_L is the window of bandpower L and Ct^ab = C^ab + N^ab. The signal C^ab comes
from forecast_3x2pt.py at the fiducial point, for all 55 pairs of the 10 tracers. The
noise N^ab is zero for a != b; for a = b it is 1/n_a for lens bins and
SIGMA_E^2/n_a for source bins, with the LSST DESC SRD Y1 number densities.

Usage: python generate_data.py
"""

import itertools
from pathlib import Path

import numpy as np
import sacc

from firecrown.likelihood import NamedParameters
from firecrown.likelihood._two_point import calculate_angular_cl
from firecrown.updatable import ParamsMap

from forecast_3x2pt import build_likelihood

HERE = Path(__file__).resolve().parent
INPUT_SACC = HERE / "paul_auto_forecast_fid_3x2pt_linear_sys_limber_20_log_bins.sacc"
OUTPUT_SACC = HERE / "forecast_3x2pt_srd_y1.sacc"
OUTPUT_NOISY_SACC = HERE / "forecast_3x2pt_srd_y1_noisy.sacc"
NOISE_SEED = 20261103

# Fiducial point. The cosmology is the DESC SRD fiducial one. The lens bias is the SRD
# Y1 bias 1.05/D(z) evaluated at the mean redshift of each lens bin, constant within the
# bin (alphag = alphaz = 0). mag_bias = 0 is the magnification slope s = 0, which keeps
# the magnification term (CCL's term is proportional to 5s - 2). ia_bias = 1 is a
# constant intrinsic alignment amplitude. Magnification and ia_bias are not SRD values.
LENS_BIAS = [1.23885511, 1.3781005, 1.52472019, 1.67665752, 1.83243479]
FIDUCIAL = {
    "Omega_c": 0.2664315,
    "Omega_b": 0.0491685,
    "h": 0.6727,
    "n_s": 0.9645,
    "sigma8": 0.831,
    "Omega_k": 0.0,
    "w0": -1.0,
    "wa": 0.0,
    "Neff": 3.046,
    "m_nu": 0.0,
    "T_CMB": 2.7255,
    "ia_bias": 1.0,
    "alphaz": 0.0,
    "z_piv": 0.0,
}
for i in range(5):
    FIDUCIAL |= {
        f"lens{i}_bias": LENS_BIAS[i],
        f"lens{i}_mag_bias": 0.0,
        f"lens{i}_delta_z": 0.0,
        f"lens{i}_alphag": 0.0,
        f"lens{i}_alphaz": 0.0,
        f"lens{i}_z_piv": 0.0,
        f"src{i}_delta_z": 0.0,
    }

# LSST DESC SRD Y1 noise: 12300 deg^2, 18 lens galaxies per arcmin^2 split over the 5
# lens bins as in Augur's srd_y1_3x2.yaml, 10 source galaxies per arcmin^2 in 5
# equal-number bins, shape noise 0.26 per component.
F_SKY = 0.3
LENS_NDENS = [2.40, 3.82, 4.29, 4.05, 3.44]  # arcmin^-2
SRC_NDENS = [2.0, 2.0, 2.0, 2.0, 2.0]  # arcmin^-2
SIGMA_E = 0.26
ARCMIN2_PER_SR = (180.0 * 60.0 / np.pi) ** 2

if __name__ == "__main__":
    # Theory vector at the fiducial point, written into a copy of the input SACC.
    sacc_in = sacc.Sacc.load_fits(str(INPUT_SACC))
    likelihood, tools = build_likelihood(NamedParameters({"sacc_file": str(INPUT_SACC)}))
    params = ParamsMap(FIDUCIAL)
    tools.update(params)
    tools.prepare()
    likelihood.update(params)
    likelihood.compute_loglike(tools)
    sacc_out = likelihood.make_realization(sacc_in, add_noise=False)

    # Signal plus noise, Ct^ab, for every pair of tracers at every ell of the windows.
    # Each C^ab is the sum over the CCL tracers of the two sources, as in TwoPoint.
    sources = {}
    for statistic in likelihood.statistics:
        for source in (statistic.statistic.source0, statistic.statistic.source1):
            sources[source.sacc_tracer] = source
    noise = {f"lens{i}": 1.0 / (LENS_NDENS[i] * ARCMIN2_PER_SR) for i in range(5)}
    noise |= {f"src{i}": SIGMA_E**2 / (SRC_NDENS[i] * ARCMIN2_PER_SR) for i in range(5)}
    window = sacc_out.data[0].tags["window"]
    ells = np.asarray(window.values, dtype=np.int64)
    ells = ells[ells >= 2]
    cl_total = {}
    for a, b in itertools.combinations_with_replacement(sorted(sources), 2):
        cl = sum(
            calculate_angular_cl(
                ells,
                f"{tracer_a.field}:{tracer_b.field}",
                sources[a].get_scale(),
                sources[b].get_scale(),
                tools,
                tracer_a,
                tracer_b,
            )
            for tracer_a, tracer_b in itertools.product(
                sources[a].get_tracers(tools), sources[b].get_tracers(tools)
            )
        )
        if a == b:
            cl = cl + noise[a]
        cl_total[(a, b)] = cl
        cl_total[(b, a)] = cl
    likelihood.reset()
    tools.reset()

    # Gaussian covariance of the bandpowers, in the data order of sacc_out. Two
    # bandpowers are correlated only if their windows overlap.
    n_data = len(sacc_out.data)
    weights = np.zeros((n_data, ells.size))
    for k, point in enumerate(sacc_out.data):
        window_k = point.tags["window"]
        assert np.array_equal(np.asarray(window_k.values)[-ells.size :], ells)
        weights[k] = window_k.weight[-ells.size :, point.tags["window_ind"]]
    variance_per_mode = 1.0 / ((2.0 * ells + 1.0) * F_SKY)
    covariance = np.zeros((n_data, n_data))
    for k in range(n_data):
        a, b = sacc_out.data[k].tracers
        for m in range(k, n_data):
            overlap = weights[k] * weights[m]
            if not overlap.any():
                continue
            c, d = sacc_out.data[m].tracers
            products = cl_total[(a, c)] * cl_total[(b, d)] + cl_total[(a, d)] * cl_total[(b, c)]
            covariance[k, m] = np.sum(overlap * products * variance_per_mode)
            covariance[m, k] = covariance[k, m]
    sacc_out.add_covariance(covariance, overwrite=True)

    # Record how the file was made, then write it.
    sacc_out.metadata["generator"] = "generate_data.py, likelihood from forecast_3x2pt.py"
    sacc_out.metadata["model"] = "pure CCL, Eisenstein-Hu and halofit, sigma8 amplitude"
    sacc_out.metadata["mean"] = "theory vector at the fiducial point, no noise"
    sacc_out.metadata["covariance"] = "Gaussian, f_sky, LSST DESC SRD Y1 noise"
    sacc_out.metadata["f_sky"] = F_SKY
    sacc_out.metadata["sigma_e"] = SIGMA_E
    for i in range(5):
        sacc_out.metadata[f"n_lens{i}_arcmin2"] = LENS_NDENS[i]
        sacc_out.metadata[f"n_src{i}_arcmin2"] = SRC_NDENS[i]
    for name, value in FIDUCIAL.items():
        sacc_out.metadata[f"fiducial_{name}"] = value
    sacc_out.save_fits(str(OUTPUT_SACC), overwrite=True)

    # Noisy data vector: the likelihood reads OUTPUT_SACC, so make_realization draws
    # the noise from the covariance computed above.
    likelihood, tools = build_likelihood(NamedParameters({"sacc_file": str(OUTPUT_SACC)}))
    tools.update(params)
    tools.prepare()
    likelihood.update(params)
    likelihood.compute_loglike(tools)
    np.random.seed(NOISE_SEED)
    sacc_noisy = likelihood.make_realization(sacc_out, add_noise=True)
    likelihood.reset()
    tools.reset()
    sacc_noisy.metadata["mean"] = "theory vector at the fiducial point plus Gaussian noise"
    sacc_noisy.metadata["noise_seed"] = NOISE_SEED
    sacc_noisy.save_fits(str(OUTPUT_NOISY_SACC), overwrite=True)
