/*
 .----------------------------------------------------------.
 |  Optimized StrayLight implementation in C++/OpenMP/CUDA  |
 |      Task 2.2 of "TASTE Maintenance and Evolutions"      |
 |           P12001 - PFL-PTE/JK/ats/412.2012i              |
 |                                                          |
 |  Contact Point:                                          |
 |           Thanassis Tsiodras, Dr.-Ing.                   |
 |           ttsiodras@semantix.gr / ttsiodras@gmail.com    |
 |           NeuroPublic S.A.                               |
 |                                                          |
 |   Licensed under the GNU Lesser General Public License,  |
 |   details here:  http://www.gnu.org/licenses/lgpl.html   |
 `----------------------------------------------------------'

*/

#include "particulate_contamination_harvey_psf.h"

//+
// NAME:
// MODIFIED_HARVEY_BRDF
//
// PURPOSE:
// Return the BRDF following a variant of the Harvey model 
//
// CATEGORY:
// Optics
//
// CALLING SEQUENCE:
// 
//
// INPUTS:
// theta: scattering angle in radian with respect to the normal of the mirror
// theta0: incidence angle in radian of main ray wrt to the normal to the surface of the mirror
// b: model parameter
// s: model parameter
// l: model parameter
//
// 
//
// OPTIONAL INPUT PARAMETERS:
// 
// OUTPUTS:
// 
// COMMON BLOCKS:
// None.
//
// SIDE EFFECTS:
// None.
//
// RESTRICTIONS:
// None.
//
// PROCEDURE:
// Originally developed by Harvey. 
//
// MODIFICATION HISTORY:
// Written mbouvet@esa.int , 2011 - AUGUST 15


fp modified_harvey_brdf(fp theta, fp theta0, fp b, fp s, fp l)
{
	return b * pow((1. + pow((sin(theta) - sin(theta0)) / l, 2)), s / 2.);
}

fp harvey_brdf(int channel, fp detector_size, int phase, fp element)
{
	// Following Peterson et al. we have ...
	// Computation of the Harvey brdf
	const fp radiusTimesNaDivMa = element * straylight::constansts::na / straylight::constansts::mirror_aperture[phase];
	const fp brdf_M =
			modified_harvey_brdf(radiusTimesNaDivMa, straylight::constansts::theta0, straylight::constansts::b_particulate_1[channel - 1],
								 straylight::constansts::s_particulate_1[channel - 1],
								 straylight::constansts::l_particulate_1[channel - 1]) +
			modified_harvey_brdf(radiusTimesNaDivMa, straylight::constansts::theta0, straylight::constansts::b_particulate_2[channel - 1],
								 straylight::constansts::s_particulate_2[channel - 1],
								 straylight::constansts::l_particulate_2[channel - 1]);
	// Irradiance in focal plane
	const fp irr_distrib_focal_M =
			sc::E_ent * M_PI * pow(straylight::constansts::mirror_aperture[0], 2) * brdf_M * straylight::constansts::na2 *
			pow(1. / straylight::constansts::mirror_aperture[phase], 2);

	// Computation of the radiant power in the focal plane at each pixel
	const fp power_focal_M = irr_distrib_focal_M * pow(detector_size, 2);

	// Computation of the normalised power distribution in focal plane
	return power_focal_M / (sc::E_ent * M_PI * pow(straylight::constansts::mirror_aperture[0], 2));
}
//+
// NAME:
// PARTICULATE_CONTAMINATION_HARVEY_PSF
//
// PURPOSE:
// Return the PSF induced by the particulate contamination of the 3 mirrors in an array of user defined dimensions
// following a variant of the Harvey model. The model paramters were obtained from fitting Mie calculation with a
// 1000 ppm contamination
//
// CATEGORY:
// Optics
//
// CALLING SEQUENCE:
// psf=particulate_CONTAMINATION_HARVEY_PSF(radius, channel)
//
// INPUTS:
// - radius: an array containing the radius to the Gaussian image in the focal plane in mm
// - channel: [1, 4]
//
// OPTIONAL INPUT PARAMETERS:
//
// OUTPUTS:
// - an array of same size than radius containing the PSF with both the direct component (1-TIS_M1-TIS_M2-TIS_M3) and the
//   diffuse component
//
// COMMON BLOCKS:
// None.
//
// SIDE EFFECTS:
// None.
//
// RESTRICTIONS:
// None.
//
// PROCEDURE:
// Model originally developed by Harvey
//
//

// MODIFICATION HISTORY:
// Written mbouvet@esa.int , 2011 - JULY 25

void particulate_contamination_harvey_psf(Matrix2D& psf_mirror_dust, const Matrix2D& radius, int channel)
{
	fp detector_size = 0.;
	if (channel <= 3) {
		detector_size = sc::VNIR_detector_size;
	}
	else if (channel == 4) {
		detector_size = sc::SWIR_detector_size;
	}



	//;;;;;;;;
	// Compute the TIS
	//;;;;;;;;;

	fp TIS_dust = 0.;
	for (int i = 0; i < 1000; ++i) {
		const fp theta_rd = float(i) / 1000. * M_PI / 2. + 1e-6;
		const fp brdf_particulate =
				modified_harvey_brdf(theta_rd, sc::theta0, sc::b_particulate_1[channel - 1], sc::s_particulate_1[channel - 1],
									 sc::l_particulate_1[channel - 1]) +
				modified_harvey_brdf(theta_rd, sc::theta0, sc::b_particulate_2[channel - 1], sc::s_particulate_2[channel - 1],
									 sc::l_particulate_2[channel - 1]);
		TIS_dust += 2. * M_PI * brdf_particulate * cos(theta_rd) * sin(theta_rd) * 1. / 1000. * M_PI / 2.;
	}

	Matrix2D mirrorDust1(psf_mirror_dust.width, psf_mirror_dust.height);
	Matrix2D mirrorDust2(psf_mirror_dust.width, psf_mirror_dust.height);
	Matrix2D mirrorDust3(psf_mirror_dust.width, psf_mirror_dust.height);

	std::transform(std::begin(radius.m_data), std::end(radius.m_data), std::begin(mirrorDust1.m_data),
				   [phase {0}, channel, detector_size](auto& element) {
					   return harvey_brdf(channel, detector_size, phase, element);
				   });
	std::transform(std::begin(radius.m_data), std::end(radius.m_data), std::begin(mirrorDust2.m_data),
				   [phase {1}, channel, detector_size](auto& element) {
					   return harvey_brdf(channel, detector_size, phase, element);
				   });
	std::transform(std::begin(radius.m_data), std::end(radius.m_data), std::begin(mirrorDust3.m_data),
				   [phase {2}, channel, detector_size](auto& element) {
					   return harvey_brdf(channel, detector_size, phase, element);
				   });

	// Add up all mirror contibutions due to surface roughness
	psf_mirror_dust.m_data = mirrorDust1.m_data + mirrorDust2.m_data + mirrorDust3.m_data;

	// Add up the direct part of the PSF
	auto minElement = std::min_element(std::begin(radius.m_data), std::end(radius.m_data));
	auto distance = std::distance(std::begin(radius.m_data), minElement);
	fp dummy_var = *minElement;
	auto jmin = distance % radius.width;
	auto imin = (distance - jmin) / radius.width;

	if (dummy_var == 0.) {
		psf_mirror_dust(imin, jmin) = pow(1. - TIS_dust, 3);
	}
}
