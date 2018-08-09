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

#include "harvey_brdf.h"

using namespace std;

fp harvey_brdf(fp theta, fp theta0, fp b, fp s, fp l, fp m, fp n)
{
	const fp g_harvey = 0.5 * (cos(theta) + cos(theta0));
	const fp f_harvey = sqrt(1.0 + pow(((sin(theta) - sin(theta0)) / l / pow(g_harvey, n)), 2));

	return b * pow(f_harvey, s) / pow(g_harvey, m);
}

void harvey_brdf(Matrix2D& result, const Matrix2D& thetaArg, fp theta0, fp b, fp s, fp l, fp m, fp n)
{
	thetaArg.clone(result);
	std::transform(std::begin(thetaArg.m_data), std::end(thetaArg.m_data), std::begin(result.m_data), [&](const auto& theta) {
		return harvey_brdf(theta, theta0, b, s, l, m, n);
	});
}
