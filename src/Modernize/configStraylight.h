#pragma once

#include "config.h"
#include <cstddef>
#include <complex>


// What precision will we use -
// single or double?
#ifdef USE_DOUBLEPRECISION
using fp = double;
using FFTWRealType= fp;
#define FFTW_PREFIX(x) fftw_ ## x
#else
using fp = float;
using FFTWRealType= fp;
#define FFTW_PREFIX(x) fftwf_ ## x
#endif

namespace straylight::constansts {
// These are the sizes computed from the IDL prototype
// They must not change across subsequent
// forward_model_single_psf_dual_resolution() calls!
constexpr size_t IMGHEIGHT = 813;
// VITO may ask us to change to 5200 pixels per scanline
constexpr size_t IMGWIDTH = 5271;
constexpr size_t KERNEL_SIZE = 11;

//# Define the surface roughness and contamination
constexpr fp surf_rough_M1 = 4e-9;
constexpr fp surf_rough_M2 = 2e-9;
constexpr fp surf_rough_M3 = 4e-9;
constexpr fp dust_ppm = 1000.;

//# pupil stop is defined later on in loop

//# SWIR ghost
constexpr int swir_ghost = 1;

constexpr fp focal = 110e-3;

//# Field of view of the instrument in degrees
constexpr fp alpha_AC_max = 17.3; // 0.25*2*0.6771;17.3 ; Max across track viewing angle in degrees
constexpr fp alpha_AL_max = 2.75; // 0.25*2*0.6771;2.75 ; Max along track viewing angle in degrees

//# Size of detectors
constexpr fp VNIR_detector_size = 13e-6; //# in meter
constexpr fp SWIR_detector_size = 25e-6; //# in meter

// Improvising a bit, to allow working in 3 modes:
constexpr int i_band = 1;

//# Pupil stop option
constexpr fp pupil_stop = 5e-10;
//if (i_band == 4) { pupil_stop=2e-9; }
//if (i_band != 4) { pupil_stop=5e-10; }

//# Derive the number of virtual detectors
const fp AL_nb_detectors = 2. * tan(alpha_AL_max * M_PI / 180.) * focal / VNIR_detector_size;
constexpr fp AC_nb_detectors = IMGWIDTH;
//if (i_band <= 3) {
//AL_nb_detectors = 2.*tan(alpha_AL_max*M_PI/180.)*focal/VNIR_detector_size;
//AC_nb_detectors = IMGWIDTH;
//}
//
//if (i_band == 4) {
//AL_nb_detectors = 2.*tan(alpha_AL_max*M_PI/180.)*focal/SWIR_detector_size;
//AC_nb_detectors = 2.*tan(alpha_AC_max*M_PI/180.)*focal/SWIR_detector_size;
//}
//numerical aperture (F#7)
constexpr fp Fnum = 7.0;
constexpr fp na = 1.0 / (2.0 * Fnum);
constexpr fp na2 = pow(na, 2);

constexpr fp wl[] = {0.45 * 1e-6, 0.65 * 1e-6, 0.825 * 1e-6, 1.6 * 1e-6};
// !!!!
//size of the instrument aperture at each mirror in meter: from Taracola code => I don't understand where these values come from
// !!!!
constexpr fp mirror_aperture[] = {9.5 * 6.0 / Fnum * 1e-3, 5.05 * 6.0 / Fnum * 1e-3, 6.85 * 6.0 / Fnum * 1e-3};

// Define the parameters for the BRDF of the mirror for each spectral band for a surface roughness of 6 nm
// These are values given to Taracola by CSL

// Blue band
constexpr fp b_6nm_blue = 82.0;
constexpr fp s_6nm_blue = -2.4;
constexpr fp l_6nm_blue = 0.005;
constexpr fp m_6nm_blue = 0.8;
constexpr fp n_6nm_blue = 0.8;

// Red band
constexpr fp b_6nm_red = 82.0;
constexpr fp s_6nm_red = -2.0;
constexpr fp l_6nm_red = 0.002;
constexpr fp m_6nm_red = 0.0;
constexpr fp n_6nm_red = 0.0;

// NIR band
constexpr fp b_6nm_nir = 30.0;
constexpr fp s_6nm_nir = -2.4;
constexpr fp l_6nm_nir = 0.0045;
constexpr fp m_6nm_nir = 0.0;
constexpr fp n_6nm_nir = 0.0;

// SWIR band
constexpr fp b_6nm_swir = 5.0;
constexpr fp s_6nm_swir = -2.5;
constexpr fp l_6nm_swir = 0.006;
constexpr fp m_6nm_swir = 0.0;
constexpr fp n_6nm_swir = 0.0;

// Put all previous values in a array for all bands
constexpr fp b_6nm[] = {b_6nm_blue, b_6nm_red, b_6nm_nir, b_6nm_swir};
constexpr fp s_6nm[] = {s_6nm_blue, s_6nm_red, s_6nm_nir, s_6nm_swir};
constexpr fp l_6nm[] = {l_6nm_blue, l_6nm_red, l_6nm_nir, l_6nm_swir};
constexpr fp m_6nm[] = {m_6nm_blue, m_6nm_red, m_6nm_nir, m_6nm_swir};
constexpr fp n_6nm[] = {n_6nm_blue, n_6nm_red, n_6nm_nir, n_6nm_swir};

// Entrance irradiance - This values has no impact on the PSF
constexpr fp E_ent = 1.; // in W.m-2

constexpr fp i_angles_center_FOV[] = {17.92, 23.03, 11.57}; // in degrees
constexpr fp i_angles_edge_FOV[] = {32.2, 42.86, 21.23};

constexpr fp delta_incidence_angle[] = {i_angles_edge_FOV[0] - i_angles_center_FOV[0], i_angles_edge_FOV[1] - i_angles_center_FOV[1],
										i_angles_edge_FOV[2] - i_angles_center_FOV[2]};
constexpr int low_to_high_spatial_resolution = 11;//  MUST BE AN ODD NUMBER!!!
}
namespace sc = straylight::constansts;