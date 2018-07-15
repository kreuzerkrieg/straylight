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

#include "harvey_psf.h"
#include "harvey_brdf.h"
#include "tis_surface_scattering_harvey.h"

//+
// NAME:
// HARVEY_PSF
//
// PURPOSE:
// Return the PSF induced by the surface roughness of the 3 mirrors in an array of user defined dimensions 
// following a variant of the Harvey model for diamond turned aluminium samples
// Surface roughness and wavelength are also input to the PSF model

// CATEGORY:
// Optics
//
// CALLING SEQUENCE:
// psf=HARVEY_PSF(radius, channel, surface_roughness_mirror1, surface_roughness_mirror2 ,surface_roughness_mirror3)
//
// INPUTS:
// - radius: an array containing the radius to the Gaussian image in the focal plane in mm
// - channel: [1, 4]
// - surface_roughness_M1 to M3: surface roughness in meter of respectively M1 to M3
//
// OPTIONAL INPUT PARAMETERS:
// - KEYWORD: INCIDENCE_ANGLE_CENTER_FOV
//            If set, this keyword set the incidence angles used for the PSF computation on each mirror to the values for the 
//            chief ray corresponding to a source in the center of the FOV
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
// Model mentioned and described in Taccola in a TN on Proba-V straylight. Largely inspired from Taracola code. 
// PSF model originally developed by Harvey  
// !!! give ref
//

// MODIFICATION HISTORY:
// Written mbouvet@esa.int , 2011 - JULY 25

using namespace std;

void PointSpreadFunction::harvey(Matrix2D& psf, const Matrix2D& radius, int channel, bool INCIDENCE_ANGLE_CENTER_FOV)
{
	//    INCIDENCE_ANGLE_CENTER_FOV=incidence_angle_center_FOV

	// Altitude
	//h=820.0

	//focal lenght
	//F=110.0

	//optical trasmission
	//const fp optics_transmission=1.0;

	fp detector_size = 0.;
	if (channel <= 3) {
		detector_size = sc::VNIR_detector_size;
	}
	if (channel == 4) {
		detector_size = sc::SWIR_detector_size;
	}

	// 
	// Total integrated scattering for 6 nm surface roughness mirrors assuming a 0 degree incidence angle
	// 
	fp TIS_6nm[4];
	TIS_6nm[0] = tis_surface_scattering_harvey(1, 0., 6e-9);
	TIS_6nm[1] = tis_surface_scattering_harvey(2, 0., 6e-9);
	TIS_6nm[2] = tis_surface_scattering_harvey(3, 0., 6e-9);
	TIS_6nm[3] = tis_surface_scattering_harvey(4, 0., 6e-9);

	// Compute the TIS for the input surface roughness

	constexpr fp surface_roughness[] = {sc::surf_rough_M1, sc::surf_rough_M2, sc::surf_rough_M3};

	// Compute the b parameter corresponding to the surface roughness
	fp b_surface_roughness[3];
	for (int i = 0; i < 3; i++) {
		const fp tmpTIS_surface_roughness = pow(4. * M_PI * surface_roughness[i] / sc::wl[channel - 1], 2);
		b_surface_roughness[i] = sc::b_6nm[channel - 1] * tmpTIS_surface_roughness / TIS_6nm[channel - 1];
	}

	//;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	// Compute the direct part of the PSF
	//;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	// When the keyword INCIDENCE_ANGLE_CENTER_FOV is set, we use the center of the FOV and the corresponding incidence
	// angles to compute the TIS

	fp theta0[3];

	if (INCIDENCE_ANGLE_CENTER_FOV) {
		theta0[0] = sc::i_angles_center_FOV[0] * M_PI / 180.;
		theta0[1] = sc::i_angles_center_FOV[1] * M_PI / 180.;
		theta0[2] = sc::i_angles_center_FOV[2] * M_PI / 180.;
	}
	else {
		theta0[0] = 0.;
		theta0[1] = 0.;
		theta0[2] = 0.;
	}

	// Compute the TIS
	fp TIS_surface_roughness[] = {tis_surface_scattering_harvey(channel, theta0[0], surface_roughness[0]),
								  tis_surface_scattering_harvey(channel, theta0[1], surface_roughness[1]),
								  tis_surface_scattering_harvey(channel, theta0[2], surface_roughness[2])};

	//;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	// Compute the diffuse part of the PSF
	//;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	// Following Peterson et al. 
	//;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


	// // Computation of the Harvey brdf
	//     Matrix2D brdf_M1;
	//     Matrix2D radiusTimesNaDivMaPlusTheta0;
	//     for(unsigned i=0; i<radius.size(); i++)
	// 	for(unsigned j=0; j<radius[0].size(); j++)
	// 	    radiusTimesNaDivMaPlusTheta0(i,j) = radius(i,j)*na/mirror_aperture[0]+theta0[0];
	//     harvey_brdf( brdf_M1, radiusTimesNaDivMaPlusTheta0, theta0[0], b_surface_roughness[0],
	// 	s_6nm[channel-1], l_6nm[channel-1], m_6nm[channel-1], n_6nm[channel-1]);
	// 
	// // Irradiance in focal plane
	//     Matrix2D psf_M1(brdf_M1.size());
	//     for(unsigned i=0; i<brdf_M1.size(); i++) {
	// 	psf_M1[i].resize(brdf_M1[0].size());
	// 	for(unsigned j=0; j<brdf_M1[0].size(); j++) {
	// 	    fp irr_distrib_focal_M1 = 
	// 		E_ent*M_PI*pow(mirror_aperture[0],2)*brdf_M1[i][j]*pow(na,2)*pow(1./mirror_aperture[0],2);
	// // Computation of the radiant power in the focal plane at each pixel
	// 	    fp power_focal_M1 = irr_distrib_focal_M1*detector_size*detector_size;
	// // Computation of the normalised power distribution in focal plane 
	// 	    fp norm_power_distrib_focal_M1 = power_focal_M1/(E_ent*M_PI*pow(mirror_aperture[0],2));
	// 	    psf_M1[i][j] = norm_power_distrib_focal_M1;
	// 	}
	//     }

	// Computation of the Harvey brdf
	psf.reset();
	Matrix2D radiusTimesNaDivMaPlusTheta0(radius.height, radius.width);
	for (int phase = 0; phase < 3; phase++) {
		radius.clone(radiusTimesNaDivMaPlusTheta0);
		radiusTimesNaDivMaPlusTheta0.m_data *= sc::na;
		radiusTimesNaDivMaPlusTheta0.m_data /= sc::mirror_aperture[phase];
		radiusTimesNaDivMaPlusTheta0.m_data += theta0[phase];

		Matrix2D brdf_M(radius.height, radius.width);
		harvey_brdf(brdf_M, radiusTimesNaDivMaPlusTheta0, theta0[phase], b_surface_roughness[phase], sc::s_6nm[channel - 1],
					sc::l_6nm[channel - 1], sc::m_6nm[channel - 1], sc::n_6nm[channel - 1]);

		// Irradiance in focal plane
		constexpr fp common = sc::E_ent * M_PI * pow(sc::mirror_aperture[0], 2);
		for (int i = 0; i < brdf_M.height; ++i) {
			for (int j = 0; j < brdf_M.width; ++j) {
				const fp irr_distrib_focal_M = common * brdf_M(i, j) * sc::na2 * pow(1. / sc::mirror_aperture[phase], 2);
				// Computation of the radiant power in the focal plane at each pixel
				const fp power_focal_M = irr_distrib_focal_M * detector_size * detector_size;
				// Computation of the normalised power distribution in focal plane 
				const fp norm_power_distrib_focal_M = power_focal_M / common;
				// Add up all mirror contibutions due to surface roughness
				psf(i, j) += norm_power_distrib_focal_M;
			}
		}
	}

	// Add up the direct part of the PSF
	auto minElement = std::min_element(std::begin(radius.m_data), std::end(radius.m_data));
	auto distance = std::distance(std::begin(radius.m_data), minElement);
	fp dummy_var = *minElement;
	auto jmin = distance % radius.width;
	auto imin = (distance - jmin) / radius.width;
	//	for (int i = 0; i < radius.height; ++i) {
	//		for (int j = 0; j < radius.width; ++j) {
	//			if (radius(i, j) < dummy_var) {
	//				imin = i;
	//		 		jmin = j;
	//				dummy_var = radius(i, j);
	//			}
	//		}
	//	}

	if (dummy_var == 0.) {
		psf(imin, jmin) = (1. - TIS_surface_roughness[0]) * (1. - TIS_surface_roughness[1]) * (1. - TIS_surface_roughness[2]);
	}
}
