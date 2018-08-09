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

#pragma once

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

#include "configStraylight.h"
#include "Matrices.h"

fp tis_surface_scattering_harvey(int channel, fp theta0, fp s_rough_mirror, bool integration_brdf);

constexpr fp tis_surface_scattering_harvey(int channel, fp theta0, fp s_rough_mirror)
{
	// Compute the TIS for the input surface roughness
	// !!!!! this expression does not result from the integration of the BRDF
	// !!!!! The modeling of the TIS variations with the incidence angle deserves further consolidation
	return pow((4. * M_PI * s_rough_mirror / sc::wl[channel - 1] * cos(theta0)), 2);
}

void tis_surface_scattering_harvey(fp* TIS, int channel, const fp* const theta0, int thetas, fp s_rough_mirror);
