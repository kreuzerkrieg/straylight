#pragma once

#include "config.h"

// These are the sizes computed from the IDL prototype
// They must not change across subsequent
// forward_model_single_psf_dual_resolution() calls!
constexpr size_t IMGHEIGHT = 813
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

