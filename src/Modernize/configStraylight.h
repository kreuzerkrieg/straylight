#pragma once

#include "config.h"

// These are the sizes computed from the IDL prototype
// They must not change across subsequent
// forward_model_single_psf_dual_resolution() calls!
constexpr size_t IMGHEIGHT = 813;
// VITO may ask us to change to 5200 pixels per scanline
constexpr size_t IMGWIDTH = 5271;
constexpr size_t KERNEL_SIZE = 11;

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
constexpr fp AL_nb_detectors = 2. * tan(alpha_AL_max * M_PI / 180.) * focal / VNIR_detector_size;
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
