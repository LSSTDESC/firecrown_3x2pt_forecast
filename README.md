# 3x2pt forecast with Firecrown

A 3x2pt analysis of a simulated LSST Y1-like data vector, set up for CosmoSIS, Cobaya
and NumCosmo with Firecrown. The data vector has 27 harmonic-space spectra from 5 lens
bins and 5 source bins: 5 galaxy clustering autospectra, 7 galaxy-galaxy lensing
spectra and 15 cosmic shear spectra, each in 20 bandpowers over 20 ≤ ℓ ≤ 2000.

## Requirements

Firecrown 1.16.1 with CosmoSIS, Cobaya and NumCosmo, for example:
```bash
conda create --name firecrown_user -c conda-forge firecrown-validated
conda activate firecrown_user
python -m pip install cobaya isitgr
```
The CosmoSIS Standard Library is not needed: the cosmology is computed by CCL inside
the likelihood.

## Files

| file | content |
|---|---|
| `forecast_3x2pt.py` | the Firecrown likelihood factory |
| `generate_data.py` | writes the two data files below |
| `generate_configs.py` | writes the CosmoSIS, Cobaya and NumCosmo configurations |
| `run_test.sh` | evaluates the likelihood once in each framework |
| `run_apes.sh` | runs the NumCosmo APES sampler under MPI (needs the `mpi_openmpi` or `mpi_mpich` build of NumCosmo) |
| `numcosmo_mpi.py` | starts the `numcosmo` command line with the Firecrown types and models registered in every MPI rank (used by `run_apes.sh`) |
| `numcosmo_catalog.py` | runs `numcosmo catalog` (`analyze`, `plot-corner`, ...) after creating the models the catalog refers to |
| `paul_auto_forecast_fid_3x2pt_linear_sys_limber_20_log_bins.sacc` | input: tracers, dN/dz and bandpower windows |
| `forecast_3x2pt_srd_y1.sacc` | noiseless data vector and covariance |
| `forecast_3x2pt_srd_y1_noisy.sacc` | the same data vector plus one realization of the noise |

## Steps

```bash
python generate_data.py                  # about 15 s
python generate_configs.py               # noisy data; add --noiseless for the noiseless file
./run_test.sh
```
`run_test.sh` prints the log-likelihood at the fiducial point. With the noisy data,
CosmoSIS and Cobaya report a log-likelihood of -269.54 and NumCosmo reports -2 times the
log-likelihood, 539.08. With
`--noiseless`, all three report 0.

`generate_configs.py` writes `cosmosis/`, `cobaya/` and `numcosmo/`. Each directory
holds its own copy of the factory and of the SACC file. The configurations evaluate the
likelihood once (CosmoSIS sampler `test`, Cobaya sampler `evaluate`,
`numcosmo run test`); change the sampler there to run a chain.

## Slides

`slides/index.qmd` is a 20-minute introduction to Firecrown and to this forecast (Quarto
RevealJS), with links to the Firecrown tutorial. It is published at
<https://lsstdesc.org/firecrown_3x2pt_forecast/>. Render it locally with
`quarto render slides/index.qmd`. On GitHub, `.github/workflows/slides.yml` renders it and
publishes it on GitHub Pages on every push to `main` that changes `slides/`; the
repository's Pages source must be set to "GitHub Actions".

## Model

`forecast_3x2pt.py` builds a Firecrown `ConstGaussian` likelihood:

* lens bin i: galaxy bias
  $b_i(z) = b_i\,[(1+z)/(1+z_{{\rm piv},i})]^{\alpha_{z,i}}\, D(z)^{\alpha_{g,i}}$ at redshift $z$,
  with $D(z)$ the linear growth factor ($D = 1$ today) (`LinearBiasSystematic`; $b_i$,
  $\alpha_{z,i}$, $z_{{\rm piv},i}$, $\alpha_{g,i}$ are `lens{i}_bias`, `lens{i}_alphaz`,
  `lens{i}_z_piv`, `lens{i}_alphag`), redshift-space distortions, magnification with
  slope s = `lens{i}_mag_bias` (`ConstantMagnificationBiasSystematic`), photo-z shift
  `lens{i}_delta_z`;
* source bin i: photo-z shift `src{i}_delta_z` and intrinsic alignment
  $A(z) = A_{\rm IA}\,[(1+z)/(1+z_{\rm piv})]^{\alpha_z}$ times $C_1\rho_{\rm crit}\Omega_m/D(z)$,
  where $C_1\rho_{\rm crit}$ is CCL's normalization of the linear alignment model
  ($\rho_{\rm crit}$ the critical density today) and $\Omega_m$ the matter density parameter
  (`LinearAlignmentSystematic`; $A_{\rm IA}$, $\alpha_z$, $z_{\rm piv}$ are `ia_bias`, `alphaz`,
  `z_piv`, shared by all source bins);
* cosmology: CCL in pure CCL mode, with the matter power spectrum P(k) from Eisenstein-Hu
  (linear) and halofit (nonlinear),
  amplitude set by `sigma8`. All three frameworks evaluate the same theory.

Spectra are computed in the Limber approximation and averaged over the top-hat
bandpower windows of the SACC file.

## Data and covariance

`generate_data.py` reads the tracers, their redshift distributions dN/dz, the list of 27 spectra and the
bandpower windows from the input SACC file. It does not use the input data vector or
covariance.

**Mean.** The theory vector of `forecast_3x2pt.py` at the fiducial point
(`generate_data.FIDUCIAL`). Ω_c, Ω_b and Ω_k are the cold dark matter, baryon and curvature
density parameters; h is the Hubble constant H₀ in units of 100 km s⁻¹ Mpc⁻¹; n_s is the scalar
spectral index; σ8 is the rms of the linear matter density contrast in spheres of radius
8 h⁻¹ Mpc today; w0 and wa set the dark energy equation of state w(a) = w0 + wa (1 - a) at scale
factor a; Σm_ν is the sum of the neutrino masses; z̄_i is the mean redshift of lens bin i.

| parameters | fiducial value | what it gives |
|---|---|---|
| Ω_c, Ω_b, h, n_s, σ8 | 0.2664315, 0.0491685, 0.6727, 0.9645, 0.831 | the LSST DESC SRD fiducial cosmology |
| w0, wa, Ω_k, Σm_ν | -1, 0, 0, 0 | flat ΛCDM with massless neutrinos |
| `lens{i}_bias` | 1.239, 1.378, 1.525, 1.677, 1.832 | 1.05/D(z̄_i), the SRD Y1 bias 1.05/D(z) evaluated at the mean redshift of lens bin i |
| `lens{i}_alphag`, `lens{i}_alphaz` | 0, 0 | bias constant in z within each lens bin |
| `lens{i}_mag_bias` | 0 | magnification included with slope s = 0: CCL's magnification term is proportional to 5s - 2, so it is not zero (it is zero for s = 0.4) |
| `ia_bias`, `alphaz`, `z_piv` | 1, 0, 0 | intrinsic alignment with constant amplitude A_IA = 1 |
| `lens{i}_delta_z`, `src{i}_delta_z` | 0 | no photo-z shifts |

The cosmology and the lens bias follow the SRD Y1 forecast; the magnification and the
intrinsic alignment amplitude are choices of this forecast, not SRD values.

**Covariance.** The data vector holds the bandpowers
$C^{ab}_B \equiv \sum_\ell w_B(\ell)\, C^{ab}_\ell$ of tracers $a$ and $b$, where $w_B$ is the
window of bandpower $B$. Their covariance $\Sigma$ is the Gaussian covariance on a fraction
$f_{\rm sky}$ of the sky,
$$
\Sigma\left[C^{ab}_B, C^{cd}_{B'}\right] = \sum_\ell \frac{w_B(\ell)\, w_{B'}(\ell)}{(2\ell+1) f_{\rm sky}} \left[\tilde C^{ac}_\ell \tilde C^{bd}_\ell + \tilde C^{ad}_\ell \tilde C^{bc}_\ell\right],
$$
where $\tilde C^{ab}_\ell \equiv C^{ab}_\ell + N^{ab}$. The noise $N^{ab}$ is zero for
$a \neq b$, $1/\bar n_a$ for a lens bin and $\sigma_e^2/\bar n_a$ for a source bin, where
$\bar n_a$ is the number density of galaxies in bin $a$ and $\sigma_e$ the shape noise per
ellipticity component. The signal $C^{ab}_\ell$ is computed for all 55 pairs of the 10 tracers
at the fiducial point. The noise follows the LSST DESC Science Requirements Document
(arXiv:1809.01669) for Y1:

| quantity | value |
|---|---|
| $f_{\rm sky}$ | 0.3 (12300 deg²) |
| lens $\bar n_a$, bins 0 to 4 | 2.40, 3.82, 4.29, 4.05, 3.44 arcmin⁻² (18 in total) |
| source $\bar n_a$ per bin | 2.0 arcmin⁻² (10 in total, equal-number bins) |
| $\sigma_e$ | 0.26 per component |

The total signal-to-noise of the noiseless data vector is 1541.

**Noise realization.** The noisy file has the mean $\mu + L \xi$, where $\mu$ is the noiseless
data vector, $\Sigma = L L^T$ is the Cholesky factorization of the covariance and $\xi$ is a
standard normal vector drawn after
`np.random.seed(20261103)` (Firecrown `make_realization`). Its chi² at the fiducial
point is 539.08 for 540 data points.

Both SACC files record these settings in their metadata.

**Tracer order of the galaxy-galaxy lensing spectra.** The 7 `galaxy_shearDensity_cl_e` spectra
store the lens tracer first (`lens0`, `src2`), as the input SACC file does; the SACC naming
convention puts the shear tracer first. Firecrown matches the spectra to their sources either
way, so the likelihood is not affected, and `firecrown sacc view --check` reports the 7 pairs. The
order is kept on purpose: the slides use it as an example of `--check`.

## Sampled parameters

`generate_configs.py` declares 8 free parameters with uniform priors: Ω_c in
[0.06, 0.46], σ8 in [0.3, 1.2], `ia_bias` in [-5, 5] and `lens{i}_bias` in [0.8, 2.5].
All other parameters are fixed at their fiducial values. The lists `CCL_PARAMETERS` and
`FIRECROWN_PARAMETERS` in `generate_configs.py` set which parameters are free and their
ranges.

## Changes to Firecrown 1.16.1 behavior

* `forecast_3x2pt.py` replaces `firecrown.likelihood._two_point.cached_angular_cl` by
  the function it wraps. The `lru_cache` hashes the CCL objects on every call, which
  takes about 30% of an evaluation, and has no hits when the cosmology changes. One
  evaluation takes about 0.7 s with 10 threads on an Intel i7-1255U.
* `generate_configs.py` replaces `firecrown.app.analysis._cosmosis.format_float`, which
  writes CosmoSIS values with 3 significant digits, by a function that writes them in
  full precision.
* `generate_configs.py` writes `cosmosis/firecrown_cosmosis.py`, which imports the
  CosmoSIS connector of the installed Firecrown, and points the `.ini` to it.
* `generate_configs.py` sets `<name>-fit: true` for the free parameters in the NumCosmo
  experiment file, where the generator writes `<name>-fit: false` for every parameter.
