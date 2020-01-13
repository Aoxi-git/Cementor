/*************************************************************************
*  2019 Janek Kozicki                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#ifndef YADE_THIN_REAL_WRAPPER_MATH_FUNCIONS_HPP
#define YADE_THIN_REAL_WRAPPER_MATH_FUNCIONS_HPP

#include <boost/cstdfloat.hpp> // Must be the first include https://www.boost.org/doc/libs/1_71_0/libs/math/doc/html/math_toolkit/rationale.html

#include <boost/config.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions.hpp>
#include <boost/math/tools/config.hpp>
#include <cmath>
#include <cstdlib>
#include <limits>

// FIXME : use correct perfect forwarding.
//         sometimes I get these errors
//             conversion of first argument to `const Real&&` would be ill-formed
//             cannot bind rvalue reference of type ‘const Real&&’ {aka ‘const ThinRealWrapper<long double>&&’} to lvalue of type ‘ThinRealWrapper<long double>’

// NOTE, these using are necessary so that a correct overload will be used regardless of what actually is Real
//		using namespace boost::multiprecision;                                                                                                         \
//		using namespace std;                                                                                                                           \

#ifdef YADE_THIN_REAL_WRAPPER_HPP
// these are inline redirections towards the correct function

#define YADE_WRAP_FUNC_1(func)                                                                                                                                 \
	inline Real func(const Real& a)                                                                                                                        \
	{                                                                                                                                                      \
		using namespace boost::multiprecision;                                                                                                         \
		using namespace std;                                                                                                                           \
		/*return func(std::forward<UnderlyingReal>(a));*/                                                                                              \
		return func(static_cast<UnderlyingReal>(a));                                                                                                   \
		/*return func(a.val);*/                                                                                                                        \
	}

#define YADE_WRAP_FUNC_2(func)                                                                                                                                 \
	inline Real func(const Real& a, const Real& b)                                                                                                         \
	{                                                                                                                                                      \
		using namespace boost::multiprecision;                                                                                                         \
		using namespace std;                                                                                                                           \
		/*return func(std::forward<UnderlyingReal>(a),std::forward<UnderlyingReal>(b));*/                                                              \
		return func(static_cast<UnderlyingReal>(a), static_cast<UnderlyingReal>(b));                                                                   \
		/*return func(a.val,b.val);*/                                                                                                                  \
	}

#else

// without ThinRealWrapper the using keywords are enough
#define YADE_WRAP_FUNC_1(func)

#define YADE_WRAP_FUNC_2(func)

#endif

// TODO: list all functions present in /usr/include/mpreal.h (from debian package libmpfrc++-dev) so we support the same set of functions as MPFR does.

// TODO: https://www.boost.org/doc/libs/1_71_0/libs/math/doc/html/math_toolkit/overview_tr1.html
// TODO: They suggest to use this -lboost_math_tr1               boost::math::acosh(x) ↔ boost::math::tr1::acosh(x)
//      ↓ …… for large scale software development where compile times are significant …… difference in performance …… as much as 20 times,
//#include <boost/math/tr1.hpp>

YADE_WRAP_FUNC_1(log)
YADE_WRAP_FUNC_1(abs)
YADE_WRAP_FUNC_1(sqrt)
YADE_WRAP_FUNC_1(sin)
YADE_WRAP_FUNC_1(cos)
YADE_WRAP_FUNC_1(acos)
YADE_WRAP_FUNC_2(pow)
YADE_WRAP_FUNC_2(atan2)

#undef YADE_WRAP_FUNC_1
#undef YADE_WRAP_FUNC_2

#endif

