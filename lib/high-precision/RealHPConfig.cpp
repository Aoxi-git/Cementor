/*************************************************************************
*  2020 Janek Kozicki                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#include <lib/high-precision/RealHPConfig.hpp>
#include <lib/high-precision/RealIO.hpp>

namespace yade {
namespace math {

	int RealHPConfig::extraStringDigits10 { 1 };

	void RealHPConfig::pyRegister()
	{
		namespace py = ::boost::python;
		py::scope cl = py::class_<RealHPConfig>(
		                       "RealHPConfig",
		                       // docstrings for static properties are forbidden in python. The solution is to put it into __doc__
		                       // https://stackoverflow.com/questions/25386370/docstrings-for-static-properties-in-boostpython
		                       "``RealHPConfig`` class provides information about RealHP<N> type.\n"
		                       ":cvar extraStringDigits10: this static variable allows o control how many extra digits to use when converting to "
		                       "decimal srings.")
		                       .add_static_property(
		                               "extraStringDigits10",
		                               py::make_getter(&RealHPConfig::extraStringDigits10, py::return_value_policy<py::return_by_value>()),
		                               py::make_setter(&RealHPConfig::extraStringDigits10, py::return_value_policy<py::return_by_value>())
		                               // python docstrings for static variables have to be written inside :cvar ………: in __doc__ of a class. See above.
		                       );
		py::def("getSupportedByEigenCgal",
		        getSupportedByEigenCgal,
		        R"""(:return: the ``tuple`` containing N from RealHP<N> precisions supported by Eigen and CGAL)""");
		py::def("getSupportedByMinieigen",
		        getSupportedByMinieigen,
		        R"""(:return: the ``tuple`` containing N from RealHP<N> precisions supported by minieigenHP)""");
		py::def("getDigits10", getDigits10, (py::arg("N")), R"""(:return: the ``int`` representing numeric_limits digits10 of RealHP<N>)""");
		py::def("getDigits2", getDigits2, (py::arg("N")), R"""(:return: the ``int`` representing numeric_limits digits (binary bits) of RealHP<N>)""");
#if (__GNUC__ < 9) // It should be checking  if (ver < 9.2.1) in fact. But cmake does the job. So here it's only to catch 'larger' mistakes.
#ifndef YADE_DISABLE_REAL_MULTI_HP
#warning "RealHP<…> won't work on this system, cmake sets YADE_DISABLE_REAL_MULTI_HP to use RealHP<1> for all precisions RealHP<N>. Also you can try -O0 flag."
// see file lib/high-precision/RealHP.hpp line: 'template <int Level> using RealHP    = Real;'
#endif
		// When using gcc older than 9.2.1 it is not possible for RealHP<N> to work. Without optimizations -O0 it can work, except for float128.
		// If YADE_DISABLE_REAL_MULTI_HP is set, then RealHP<1> is used in place of all possible precisions RealHP<N> : see file RealHP.hpp for this setting.
		// So this is for local testing only. With flag -O0 most of RealHP<…> works, except for float128 which is always segfaulting.
		py::scope().attr("isFloat128Broken") = true;
#else
		py::scope().attr("isFloat128Broken") = false;
#endif
#ifndef YADE_DISABLE_REAL_MULTI_HP
		py::scope().attr("isEnabledRealHP") = true;
#else
		py::scope().attr("isEnabledRealHP") = false;
#endif
		py::scope().attr("workaroundSlowBoostBinFloat") = int(workaroundSlowBoostBinFloat);
	}

} // namespace math
} // namespace yade

/* FIXME - put this into documentation (above) of extraStringDigits10

The extraStringDigits10 is to make sure that there are no conversion errors in the last bit.
here is a quick python example which shows the 'cutting' of last digits.

# This one demostrates that `double` used by python has just 53 bits of precision:

for a in range(128): print(str(a).rjust(3,' ')+": "+str(1+1./pow(2.,a)))

# This one shows the 'true' values:

import mpmath; mpmath.mp.dps=200;
for a in range(128): print(str(a).rjust(3,' ')+": "+str(mpmath.mpf(1)+mpmath.mpf(1)/pow(mpmath.mpf(2),mpmath.mpf(a))))

# This one shows the actual 'Real' precision used in yade. To achieve this the mth.max(…,…) are called to force the numbers
# to pass through C++, instead of letting mpmath to calculate this, so for example we can see that float128 has 113 bits.
# Also this test was used to verify the value for extraStringDigits10 as well as the formula given
# in IEEE Std 754-2019 Standard for Floating-Point Arithmetic: Pmin (bf) = 1 + ceiling( p × log10(2)), where p is the number of significant bits in bf

from yade import math as mth
for a in range(128): print(str(a).rjust(3,' ')+": "+str(mth.max(0,mth.max(0,1)+mth.max(0,1)/mth.pow(mth.max(0,2),a))))

# But also we can now check the precision directly by calling

yade.math.RealHPConfig.getDigits2(N) # for any N of RealHP<N>

# Hmm maybe turn this into an external parameter? Configurable from python? And write in help "use 1 to extract results and avoid fake sense of more precision,
# use 4 or more to have numbers which will convert exactly in both directions mpmath ↔ string ↔ C++.".
# For now it is in yade.math.RealHPConfig.extraStringDigits10
*/
