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

#include "tis_surface_scattering_harvey.h"
#include "harvey_brdf.h"
#include "utilities.h"

using namespace std;

//+
// NAME:
// TIS_SURFACE_SCATTERING_HARVEY
//
// PURPOSE:
// Return the Total Integrated Scattering for the Harvey model for diamond turned aluminium samples, by default
// this is done using a cos(theta_incident)^2 law. The computation of the TIS can however also be done on the BRDF directly, by integration. 
//
// CATEGORY:
// Optics
//
// CALLING SEQUENCE:
// 
//
// INPUTS:
// 
// - theta0: incidence angle of main ray wrt to the normal to the surface of the mirror
// - channel: 1 to 4
// - s_rough_mirror: surface roughness in nm
// 
//
// OPTIONAL INPUT PARAMETERS:
// - /INTEGRATION_BRDF : to allow for direct integration of the Harvey BRDF
// OUTPUTS:
// - scalar or vector TIS : same dimensions as theta0


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
// Model mentioned and described in Taccola in a TN on Proba-V straylight. Originally developed by Harvey. 
//
// MODIFICATION HISTORY:
// Written mbouvet@esa.int , 2011 - AUGUST 25


fp tis_surface_scattering_harvey(int channel, fp theta0, fp s_rough_mirror, bool integration_brdf)
{
	// Compute the TIS for the input surface roughness
	// !!!!! this expression does not result from the integration of the BRDF
	// !!!!! The modeling of the TIS variations with the incidence angle deserves further consolidation
	fp TIS = pow((4. * M_PI * s_rough_mirror / sc::wl[channel - 1] * cos(theta0)), 2);

	if (integration_brdf) {
		debug_printf(LVL_WARN, "Untested code, (Line %d in %s)\n", __LINE__, __FILE__);

		// Integration step
		constexpr int nb_points = 1200; // 600 is a good values to get the fwd models to run reasonably fast but not completely accurate

		// compute the TIS for a normal incident angle
		fp TIS_6nm = 0.;

		for (int i = 0; i < nb_points; ++i) {
			const fp theta_rd = (fp(i) / nb_points) * M_PI / 2.;
			TIS_6nm += 2. * M_PI * harvey_brdf(theta_rd, 0, sc::b_6nm[channel - 1], sc::s_6nm[channel - 1], sc::l_6nm[channel - 1],
											   sc::m_6nm[channel - 1], sc::n_6nm[channel - 1]) * cos(theta_rd) * sin(theta_rd) * 1. /
					   nb_points * M_PI / 2.;
		}

		// Compute the TIS for the input surface roughness
		const fp TIS_surface_roughness_norm = pow(4. * M_PI * s_rough_mirror / sc::wl[channel - 1], 2);

		// Compute the b parameter corresponding to the surface roughness
		const fp b_surface_roughness = sc::b_6nm[channel - 1] * TIS_surface_roughness_norm / TIS_6nm;

		TIS = 0.;

		for (int i_phi = 0; i_phi < nb_points; i_phi++) {
			fp phi_rd = (fp(i_phi) / nb_points) * M_PI / 2.;
			fp TISdelta = 0.;
			for (int i = 0; i < nb_points; ++i) {
				const fp theta_rd = (fp(i) / nb_points) * M_PI / 2.;
				fp tmp = (cos(theta_rd) * cos(theta0) + sin(theta_rd) * sin(theta0) * cos(phi_rd));
				if (tmp > 1.) {
					tmp = 1.;
				}
				if (tmp < -1.) {
					tmp = -1.;
				}
				fp scattering_angle_wrt_incoming = acos(tmp);

				fp hv = harvey_brdf(-scattering_angle_wrt_incoming + theta0, theta0, b_surface_roughness, sc::s_6nm[channel - 1],
									sc::l_6nm[channel - 1], sc::m_6nm[channel - 1], sc::n_6nm[channel - 1]);
				TISdelta += hv * cos(theta_rd) * sin(theta_rd) * 1. / nb_points * M_PI / 2.;
			}
			TIS += 1. / nb_points * M_PI * TISdelta;
		}
		// We have integrated only on phi from 0 to 180. The BRDF is symmetrical so we double the TIS:
		TIS *= 2.;
	}
	return TIS;
}

void tis_surface_scattering_harvey(fp* TIS, int channel, const fp* const theta0, int thetas, fp s_rough_mirror)
{
	// Compute the TIS for the input surface roughness
	// !!!!! this expression does not result from the integration of the BRDF
	// !!!!! The modeling of the TIS variations with the incidence angle deserves further consolidation
	for (int i = 0; i < thetas; ++i) {
		TIS[i] = tis_surface_scattering_harvey(channel, theta0[i], s_rough_mirror);
	}
}
/*void tis_surface_scattering_harvey(fp* TIS, int channel, Matrix1D& theta0, fp s_rough_mirror, bool integration_brdf)
{
	// Compute the TIS for the input surface roughness
	// !!!!! this expression does not result from the integration of the BRDF
	// !!!!! The modeling of the TIS variations with the incidence angle deserves further consolidation
	for (int i = 0; i < theta0.m_data.size(); ++i) {
		TIS[i] = (4. * M_PI * s_rough_mirror / sc::wl[channel - 1] * cos(theta0[i]));
		TIS[i] *= TIS[i];
	}

	if (integration_brdf) {
		debug_printf(LVL_WARN, "Untested code, (Line %d in %s)\n", __LINE__, __FILE__);

		// Integration step
		constexpr int nb_points = 1200; // 600 is a good values to get the fwd models to run reasonably fast but not completely accurate

		// compute the TIS for a normal incident angle
		fp TIS_6nm = 0.;

		for (int i = 0; i < nb_points; ++i) {
			const fp theta_rd = (fp(i) / nb_points) * M_PI / 2.;
			TIS_6nm += 2. * M_PI * harvey_brdf(theta_rd, 0, sc::b_6nm[channel - 1], sc::s_6nm[channel - 1], sc::l_6nm[channel - 1],
											   sc::m_6nm[channel - 1], sc::n_6nm[channel - 1]) * cos(theta_rd) * sin(theta_rd) * 1. /
					   nb_points * M_PI / 2.;
		}

		// Compute the TIS for the input surface roughness
		const fp TIS_surface_roughness_norm = pow(4. * M_PI * s_rough_mirror / sc::wl[channel - 1], 2);

		// Compute the b parameter corresponding to the surface roughness
		const fp b_surface_roughness = sc::b_6nm[channel - 1] * TIS_surface_roughness_norm / TIS_6nm;

		//memset(TIS, 0, sizeof(fp) * theta0.m_data.size());
		const int nb_i_angles = theta0.m_data.size();
		for (int i_angles = 0; i_angles < nb_i_angles; ++i_angles) {
			for (int i_phi = 0; i_phi < nb_points; ++i_phi) {
				fp phi_rd = (fp(i_phi) / nb_points) * M_PI / 2.;
				fp TISdelta = 0.;
				for (int i = 0; i < nb_points; ++i) {
					fp theta_rd = (fp(i) / nb_points) * M_PI / 2.;
					fp tmp = (cos(theta_rd) * cos(theta0[i_angles]) + sin(theta_rd) * sin(theta0[i_angles]) * cos(phi_rd));
					if (tmp > 1.) {
						tmp = 1.;
					}
					if (tmp < -1.) {
						tmp = -1.;
					}
					fp scattering_angle_wrt_incoming = acos(tmp);

					fp hv = harvey_brdf(-scattering_angle_wrt_incoming + theta0[i_angles], theta0[i_angles], b_surface_roughness,
										sc::s_6nm[channel - 1], sc::l_6nm[channel - 1], sc::m_6nm[channel - 1], sc::n_6nm[channel - 1]);
					TISdelta += hv * cos(theta_rd) * sin(theta_rd) * 1. / nb_points * M_PI / 2.;
				}
				TIS[i_angles] += 1. / nb_points * M_PI * TISdelta;
			}
		}
		// We have integrated only on phi from 0 to 180. The BRDF is symmetrical so we double the TIS:
		for (int i_angles = 0; i_angles < nb_i_angles; ++i_angles) {
			TIS[i_angles] *= 2.;
		}
	}
}*/
