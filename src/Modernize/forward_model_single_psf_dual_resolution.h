#pragma once

#include "Matrices.h"
#include "configStraylight.h"
#include <fftw3.h>

using m2dSpecialPSF = Matrix2D;

class FMSPSFDR
{
public:
	void forward_model_single_psf_dual_resolution(Matrix2D& output_image, const Matrix2D& input_image, int channel, bool logStages = false);

private:
	void calculatePSF(int channel);
	//	constexpr auto image_dim_x = sc::IMGWIDTH;
	//	constexpr auto image_dim_y = sc::IMGHEIGHT;

	Matrix1D correction_direct_peak {sc::IMGWIDTH};

	static constexpr auto tmp_extended_image_dim_x = 3 * sc::IMGWIDTH;
	static constexpr auto tmp_extended_image_dim_y = 3 * sc::IMGHEIGHT;
	//#
	//#  Iterate on dimensions that are odd and multiple of sc::low_to_high_spatial_resolution and larger than 3*sc::IMGWIDTH
	//#
	static constexpr auto cpt_dim_y = []() constexpr {
		size_t retVal = 0;
		do {
			retVal++;
		} while (!(((tmp_extended_image_dim_y + retVal) % 2 == 1) &&
				   ((tmp_extended_image_dim_y + retVal) % sc::low_to_high_spatial_resolution == 0)));
		return retVal;
	}();
	static constexpr auto cpt_dim_x = []() constexpr {
		size_t retVal = 0;
		do {
			retVal++;
		} while (!(((tmp_extended_image_dim_x + retVal) % 2 == 1) &&
				   ((tmp_extended_image_dim_x + retVal) % sc::low_to_high_spatial_resolution == 0)));
		return retVal;
	}();

	//#  Define extended input image which is used during the FFT
	//#  This extended image has dimensions which are twice the input_image passed to the routine
	//#  The input_image variable is place at the center of the newly created variable
	//#  The newly created variable must be of dimension at least 3 x dimensions of the original image and should have an odd number of
	//#  columns and lines and should be a multiple of sc::low_to_high_spatial_resolution
	//#
	//#  Define starting values of the dimensions of the extended image:
	static constexpr auto extended_image_dim_x = 3 * sc::IMGWIDTH + cpt_dim_x;
	static constexpr auto extended_image_dim_y = 3 * sc::IMGHEIGHT + cpt_dim_y;

	static constexpr auto lowres_y = extended_image_dim_y / sc::low_to_high_spatial_resolution;
	static constexpr auto lowres_x = extended_image_dim_x / sc::low_to_high_spatial_resolution;
	int psf_extent = 0;
	m2dSpecialPSF* psf_high_res = nullptr;
	FFTW_PREFIX(complex)* fft_psf_low_res = nullptr;
};