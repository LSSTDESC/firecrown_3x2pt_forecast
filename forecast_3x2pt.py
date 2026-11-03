"""Firecrown likelihood factory for a 3x2pt forecast with 5 lens and 5 source bins.

build_likelihood returns a Gaussian likelihood of 27 harmonic-space spectra: 5 galaxy
clustering autospectra, 7 galaxy-galaxy lensing spectra (GGL_PAIRS, source bin behind
lens bin) and 15 cosmic shear spectra.

Lens bin i has the galaxy bias of LinearBiasSystematic,
lens{i}_bias ((1 + z)/(1 + lens{i}_z_piv))^lens{i}_alphaz D(z)^lens{i}_alphag,
redshift-space distortions, magnification with slope s = lens{i}_mag_bias (CCL's
magnification term is proportional to 5s - 2, so s = 0 does not remove it) and a photo-z
shift lens{i}_delta_z. Source bin i has a photo-z
shift src{i}_delta_z and the intrinsic alignment of LinearAlignmentSystematic, whose
parameters ia_bias, alphaz and z_piv are shared by all source bins.

CCL computes the cosmology in pure CCL mode: Eisenstein-Hu linear P(k), halofit
nonlinear P(k), amplitude set by sigma8. CosmoSIS, Cobaya and NumCosmo therefore all
evaluate the same theory.

Importing this module replaces Firecrown's cached_angular_cl by the function it wraps.
Its lru_cache hashes the CCL objects on every call, which takes about 30% of one
likelihood evaluation, and it has no hits when the cosmology changes between calls.
"""

import os

import firecrown.likelihood._two_point as _two_point
import firecrown.likelihood.number_counts as nc
import firecrown.likelihood.weak_lensing as wl
from firecrown.likelihood import ConstGaussian, Likelihood, NamedParameters, TwoPoint
from firecrown.likelihood.factories import load_sacc_data
from firecrown.modeling_tools import CCLFactory, ModelingTools
from firecrown.modeling_tools._ccl_factory import (
    CCLCreationMode,
    CCLPureModeTransferFunction,
    PoweSpecAmplitudeParameter,
)

# getattr keeps the replacement correct when the module is loaded twice in one process.
_two_point.cached_angular_cl = getattr(
    _two_point.cached_angular_cl, "__wrapped__", _two_point.cached_angular_cl
)

N_LENS = 5
N_SRC = 5
GGL_PAIRS = [(0, 2), (0, 3), (0, 4), (1, 3), (1, 4), (2, 4), (3, 4)]  # (lens, source)


def build_likelihood(params: NamedParameters) -> tuple[Likelihood, ModelingTools]:
    """Return the 3x2pt Gaussian likelihood and its modeling tools.

    :param params: build parameters; sacc_file is the path of the SACC file to read
        (environment variables and ~ are expanded)
    :returns: the ConstGaussian likelihood holding the 27 spectra in the order
        clustering, galaxy-galaxy lensing, cosmic shear, and the ModelingTools
    """
    sacc_data = load_sacc_data(
        os.path.expanduser(os.path.expandvars(params.get_string("sacc_file")))
    )

    lenses = [
        nc.NumberCounts(
            sacc_tracer=f"lens{i}",
            has_rsd=True,
            systematics=[
                nc.PhotoZShift(sacc_tracer=f"lens{i}"),
                nc.ConstantMagnificationBiasSystematic(sacc_tracer=f"lens{i}"),
                nc.LinearBiasSystematic(sacc_tracer=f"lens{i}"),
            ],
        )
        for i in range(N_LENS)
    ]
    # One alignment systematic shared by all source bins: ia_bias, alphaz and z_piv
    # are global parameters.
    ia_systematic = wl.LinearAlignmentSystematic()
    sources = [
        wl.WeakLensing(
            sacc_tracer=f"src{i}",
            systematics=[wl.PhotoZShift(sacc_tracer=f"src{i}"), ia_systematic],
        )
        for i in range(N_SRC)
    ]

    statistics = [
        TwoPoint(source0=lenses[i], source1=lenses[i], sacc_data_type="galaxy_density_cl")
        for i in range(N_LENS)
    ]
    statistics += [
        TwoPoint(
            source0=lenses[i], source1=sources[j], sacc_data_type="galaxy_shearDensity_cl_e"
        )
        for i, j in GGL_PAIRS
    ]
    statistics += [
        TwoPoint(source0=sources[i], source1=sources[j], sacc_data_type="galaxy_shear_cl_ee")
        for i in range(N_SRC)
        for j in range(i, N_SRC)
    ]

    likelihood = ConstGaussian(statistics=statistics)
    likelihood.read(sacc_data)

    modeling_tools = ModelingTools(
        ccl_factory=CCLFactory(
            amplitude_parameter=PoweSpecAmplitudeParameter.SIGMA8,
            creation_mode=CCLCreationMode.PURE_CCL_MODE,
            pure_ccl_transfer_function=CCLPureModeTransferFunction.EISENSTEIN_HU,
        )
    )
    return likelihood, modeling_tools
