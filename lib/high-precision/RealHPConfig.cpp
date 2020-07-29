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

	template <template <int> class dig> int RealHPConfig::getDigits(int N)
	{
		// 5 is the largest length of RealHPLadder<…>. If more were added, and precision were not the multiplies of digits10*N
		// then the python test will quickly catch that problem. Then more cases will be needed to add to this switch.
		static_assert(
		        boost::mpl::size<::yade::math::RealHPLadder>::value <= 5,
		        "More types were added in RealHP.hpp, please adjust this switch(…) accordingly.");
		switch (N) {
			// NOTE: I could turn this into a self-unfolding macro, like in RealHPEigenCgal.hpp (there it was necessary), but maybe better to have here a simple readable code?
			case 1: return dig<1>::value();
			case 2: return dig<2>::value();
			case 3: return dig<3>::value();
			case 4: return dig<4>::value();
			case 5: return dig<5>::value();
			case 6: return dig<6>::value();
			case 7: return dig<7>::value();
			case 8: return dig<8>::value();
			case 9: return dig<9>::value();
			case 10: return dig<10>::value();
			case 20: return dig<20>::value();
			default: return dig<1>::value() * N; // this formula is used by NthLevel in lib/high-precision/RealHP.hpp
		}
	}

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

