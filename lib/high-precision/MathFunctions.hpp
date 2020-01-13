/*************************************************************************
*  2019 Janek Kozicki                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

// This file contains mathematical functions available in standard library and boost library.
//     https://en.cppreference.com/w/cpp/numeric/math
//     https://en.cppreference.com/w/cpp/numeric/special_functions
// They have to be provided here as inline redirections towards the correct implementation, depending on what precision type yade is being compiled with.
// This is the only way to make sure that ::std, ::boost::math, ::boost::multiprecision are all called correctly.

// TODO: Boost documentation recommends to link with tr1: -lboost_math_tr1 as it provides significant speedup. For example replace boost::math::acosh(x) ↔ boost::math::tr1::acosh(x)
//     https://www.boost.org/doc/libs/1_71_0/libs/math/doc/html/math_toolkit/overview_tr1.html
//#include <boost/math/tr1.hpp>


#ifndef YADE_THIN_REAL_WRAPPER_MATH_FUNCIONS_HPP
#define YADE_THIN_REAL_WRAPPER_MATH_FUNCIONS_HPP

#include <boost/config.hpp>
#include <boost/math/complex.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions.hpp>
#include <boost/math/tools/config.hpp>
#include <boost/random.hpp>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <limits>
#include <utility>

#ifndef YADE_REAL_MATH_NAMESPACE
#error "This file cannot be included alone, include Real.hpp instead"
#endif

// Macors for quick inline redirections towards the correct function from (1) standard library or (2) boost::multiprecision; depending on which one is used.
#define YADE_WRAP_FUNC_1(func)                                                                                                                                 \
	inline Real func(const Real& a)                                                                                                                        \
	{                                                                                                                                                      \
		using YADE_REAL_MATH_NAMESPACE::func;                                                                                                          \
		using ::std::func;                                                                                                                             \
		return func(static_cast<const UnderlyingReal&>(a));                                                                                            \
	}

#define YADE_WRAP_FUNC_1_RENAME(func1, func2)                                                                                                                  \
	inline Real func1(const Real& a) { return YADE_REAL_MATH_NAMESPACE::func2(static_cast<const UnderlyingReal&>(a)); }

#define YADE_WRAP_FUNC_2(func)                                                                                                                                 \
	inline Real func(const Real& a, const Real& b)                                                                                                         \
	{                                                                                                                                                      \
		return YADE_REAL_MATH_NAMESPACE::func(static_cast<const UnderlyingReal&>(a), static_cast<const UnderlyingReal&>(b));                           \
	}

#define YADE_WRAP_FUNC_2_TYPE2(func, SecondType)                                                                                                               \
	inline Real func(const Real& a, SecondType b) { return YADE_REAL_MATH_NAMESPACE::func(static_cast<const UnderlyingReal&>(a), b); }

#define YADE_WRAP_FUNC_2_TYPE2_STD_CAST(func, SecondType)                                                                                                      \
	inline Real func(const Real& a, SecondType b) { return ::std::func(static_cast<const UnderlyingReal&>(a), static_cast<const UnderlyingReal&>(b)); }

#ifdef YADE_THIN_REAL_WRAPPER_HPP
#define YADE_WRAP_FUNC_2_TYPE2_CAST(func, SecondType, CastType)                                                                                                \
	inline Real func(const Real& a, SecondType b) { return YADE_REAL_MATH_NAMESPACE::func(static_cast<const UnderlyingReal&>(a), b->operator CastType()); }
#else
#define YADE_WRAP_FUNC_2_TYPE2_CAST(func, SecondType, CastType)                                                                                                \
	inline Real func(const Real& a, SecondType b) { return YADE_REAL_MATH_NAMESPACE::func(static_cast<const UnderlyingReal&>(a), b); }
#endif

#define YADE_WRAP_FUNC_3(func)                                                                                                                                 \
	inline Real func(const Real& a, const Real& b, const Real& c)                                                                                          \
	{                                                                                                                                                      \
		return YADE_REAL_MATH_NAMESPACE::func(                                                                                                         \
		        static_cast<const UnderlyingReal&>(a), static_cast<const UnderlyingReal&>(b), static_cast<const UnderlyingReal&>(c));                  \
	}

#define YADE_WRAP_FUNC_3_TYPE3(func, ThirdType)                                                                                                                \
	inline Real func(const Real& a, const Real& b, ThirdType c)                                                                                            \
	{                                                                                                                                                      \
		return YADE_REAL_MATH_NAMESPACE::func(static_cast<const UnderlyingReal&>(a), static_cast<const UnderlyingReal&>(b), c);                        \
	}

#define YADE_WRAP_FUNC_1_COMPLEX(func)                                                                                                                         \
	inline Complex func(const Complex& a)                                                                                                                  \
	{                                                                                                                                                      \
		using YADE_REAL_MATH_NAMESPACE::func;                                                                                                          \
		using ::std::func;                                                                                                                             \
		return func(static_cast<const std::complex<UnderlyingReal>&>(a));                                                                              \
	}

#define YADE_WRAP_FUNC_1_COMPLEX_STD(func)                                                                                                                     \
	inline Complex func(const Complex& a) { return ::std::func(static_cast<const std::complex<UnderlyingReal>&>(a)); }

#define YADE_WRAP_FUNC_1_COMPLEX_TO_REAL(func)                                                                                                                 \
	inline Real func(const Complex& a)                                                                                                                     \
	{                                                                                                                                                      \
		using YADE_REAL_MATH_NAMESPACE::func;                                                                                                          \
		using ::std::func;                                                                                                                             \
		return func(static_cast<const std::complex<UnderlyingReal>&>(a));                                                                              \
	}

#define YADE_WRAP_FUNC_1_COMPLEX_TO_REAL_STD(func)                                                                                                             \
	inline Real func(const Complex& a) { return ::std::func(static_cast<const std::complex<UnderlyingReal>&>(a)); }


namespace yade {
namespace math {
	/********************************************************************************************/
	/**********************            trigonometric functions             **********************/
	/********************************************************************************************/

	YADE_WRAP_FUNC_1(sin)
	YADE_WRAP_FUNC_1(sinh)
	YADE_WRAP_FUNC_1(cos)
	YADE_WRAP_FUNC_1(cosh)
	YADE_WRAP_FUNC_1(tan)
	YADE_WRAP_FUNC_1(tanh)

	/**********************                    Complex                     *********************/
	// add more complex functions as necessary, but remember to add them in py/high-precision/_math.cpp, py/tests/testMath.py and py/tests/testMathHelper.py
	YADE_WRAP_FUNC_1_COMPLEX(sin)
	YADE_WRAP_FUNC_1_COMPLEX(sinh)
	YADE_WRAP_FUNC_1_COMPLEX(cos)
	YADE_WRAP_FUNC_1_COMPLEX(cosh)
	YADE_WRAP_FUNC_1_COMPLEX(tan)
	YADE_WRAP_FUNC_1_COMPLEX(tanh)

	/********************************************************************************************/
	/**********************        inverse trigonometric functions         **********************/
	/********************************************************************************************/

	YADE_WRAP_FUNC_1(asin)
	YADE_WRAP_FUNC_1(asinh)
	YADE_WRAP_FUNC_1(acos)
	YADE_WRAP_FUNC_1(acosh)
	YADE_WRAP_FUNC_1(atan)
	YADE_WRAP_FUNC_1(atanh)
	YADE_WRAP_FUNC_2(atan2)

	/********************************************************************************************/
	/**********************   logarithm, exponential and power functions   **********************/
	/********************************************************************************************/

	YADE_WRAP_FUNC_1(log)
	YADE_WRAP_FUNC_1(log10)
	YADE_WRAP_FUNC_1(log1p)
	YADE_WRAP_FUNC_1(log2)
	YADE_WRAP_FUNC_1(logb)
	YADE_WRAP_FUNC_1(ilogb)

	YADE_WRAP_FUNC_2_TYPE2(ldexp, int)
	YADE_WRAP_FUNC_2_TYPE2(frexp, int*) // that's original C signature of this function

	YADE_WRAP_FUNC_1(exp)
	YADE_WRAP_FUNC_1(exp2)
	YADE_WRAP_FUNC_1(expm1)

	YADE_WRAP_FUNC_2(pow)
	YADE_WRAP_FUNC_1(sqrt)
	YADE_WRAP_FUNC_1(cbrt)

	YADE_WRAP_FUNC_2(hypot)
	//YADE_WRAP_FUNC_3(hypot) // since C++17, could be very useful for us

	/**********************                    Complex                     *********************/
	// add more complex functions as necessary, but remember to add them in py/high-precision/_math.cpp, py/tests/testMath.py and py/tests/testMathHelper.py
	YADE_WRAP_FUNC_1_COMPLEX(exp)
	YADE_WRAP_FUNC_1_COMPLEX(log)

	/********************************************************************************************/
	/**********************    min, max, abs, sign, floor, ceil, etc...    **********************/
	/********************************************************************************************/

	// Both must be found by automatic lookup: the ones from ::std and the ones that accept non-double Real types.
	using ::std::abs;
	using ::std::fabs;
	using ::std::max;
	using ::std::min;
#if (YADE_REAL_BIT != 64)
	YADE_WRAP_FUNC_2_TYPE2_STD_CAST(min, const double&)
	YADE_WRAP_FUNC_2_TYPE2_STD_CAST(max, const double&)
	YADE_WRAP_FUNC_2_TYPE2_STD_CAST(min, const Real&)
	YADE_WRAP_FUNC_2_TYPE2_STD_CAST(max, const Real&)
#endif
#if (YADE_REAL_BIT > 64)
	YADE_WRAP_FUNC_1(abs)
	YADE_WRAP_FUNC_1_RENAME(fabs, abs)
#endif
	template <typename T> int sgn(T val) { return (T(0) < val) - (val < T(0)); }
	template <typename T> int sign(T val) { return (T(0) < val) - (val < T(0)); }
	YADE_WRAP_FUNC_1(floor)
	YADE_WRAP_FUNC_1(ceil)
	YADE_WRAP_FUNC_1(round)
	YADE_WRAP_FUNC_1(rint)
	YADE_WRAP_FUNC_1(trunc)

#ifndef YADE_IGNORE_IEEE_INFINITY_NAN
	YADE_WRAP_FUNC_1(isnan)
	YADE_WRAP_FUNC_1(isinf)
	YADE_WRAP_FUNC_1(isfinite)
#endif

	/**********************                    Complex                     *********************/
	// add more complex functions as necessary, but remember to add them in py/high-precision/_math.cpp and py/tests/testMath.py
	YADE_WRAP_FUNC_1_COMPLEX_STD(conj)
	YADE_WRAP_FUNC_1_COMPLEX_TO_REAL(abs)
	YADE_WRAP_FUNC_1_COMPLEX_TO_REAL_STD(real)
	YADE_WRAP_FUNC_1_COMPLEX_TO_REAL_STD(imag)

	/********************************************************************************************/
	/**********************        integer division and remainder          **********************/
	/********************************************************************************************/

	YADE_WRAP_FUNC_2(fmod)
	YADE_WRAP_FUNC_2(remainder)

	YADE_WRAP_FUNC_2_TYPE2_CAST(modf, Real*, UnderlyingReal*)

	YADE_WRAP_FUNC_3(fma)
	YADE_WRAP_FUNC_3_TYPE3(remquo, int*)

	/********************************************************************************************/
	/**********************         special mathematical functions         **********************/
	/********************************************************************************************/

	YADE_WRAP_FUNC_1(erf)
	YADE_WRAP_FUNC_1(erfc)
	YADE_WRAP_FUNC_1(lgamma)

// These will be available in C++17, we could use the ones from boost, if they become necessary.
//YADE_WRAP_FUNC_1(riemann_zeta)
//YADE_WRAP_FUNC_2(beta)
//YADE_WRAP_FUNC_2(cyl_bessel_i)
//YADE_WRAP_FUNC_2(cyl_bessel_j)
//YADE_WRAP_FUNC_2(cyl_bessel_k)
//YADE_WRAP_FUNC_2_TYPE1(sph_bessel, unsigned)


// workaround broken tgamma for boost::float128
#if (YADE_REAL_BIT <= 128) and (YADE_REAL_BIT > 80)
	static_assert(std::is_same<UnderlyingReal, boost::multiprecision::float128>::value, "Incorrect type, please file a bug report.");
	inline Real tgamma(const Real& a)
	{
		if (a >= 0) {
			return YADE_REAL_MATH_NAMESPACE::tgamma(static_cast<UnderlyingReal>(a));
		} else {
			return abs(YADE_REAL_MATH_NAMESPACE::tgamma(static_cast<UnderlyingReal>(a)))
			        * ((static_cast<unsigned long long>(floor(abs(a))) % 2 == 0) ? -1 : 1);
		}
	}
#else
	YADE_WRAP_FUNC_1(tgamma)
#endif

	/********************************************************************************************/
	/**********************        extract C-array from std::vector        **********************/
	/********************************************************************************************/

// Some old C library functions need pointer to C-array, this is for compatibility between ThinRealWrapper and UnderlyingReal
#ifdef YADE_THIN_REAL_WRAPPER_HPP
	static_assert(sizeof(Real) == sizeof(UnderlyingReal), "This compiler introduced padding. This breaks binary compatibility");
	static_assert(sizeof(Complex) == sizeof(std::complex<UnderlyingReal>), "This compiler introduced padding, which breaks binary compatibility");

	static inline const UnderlyingReal* constVectorData(const std::vector<Real>& v) { return v.data()->operator const UnderlyingReal*(); }
	static inline UnderlyingReal*       vectorData(std::vector<Real>& v) { return v.data()->operator UnderlyingReal*(); }
#else
	static inline const UnderlyingReal* constVectorData(const std::vector<Real>& v) { return v.data(); }
	static inline UnderlyingReal*       vectorData(std::vector<Real>& v) { return v.data(); }
#endif

	/********************************************************************************************/
	/**********************                     random                     **********************/
	/********************************************************************************************/

	// These random functions are necessary for Eigen library to for example write in python: Vector3.Random()
	// generate random number [0,1)
	static inline Real random01()
	{
		static ::boost::random::mt19937 gen;
		return ::boost::random::generate_canonical<::yade::math::Real, std::numeric_limits<::yade::math::Real>::digits>(gen);
	}

	static inline Real unitRandom() { return random01(); }
	static inline Real random() { return random01() * 2 - 1; }

}
// do we want such alias?
//namespace m = ::yade::math;
//namespace mth = ::yade::math;

// Use those which argument dependent lookup (ADL) can't find, because sometimes they are used with non-Real type, like int or float.
// This saves writing `math::` before the calls. If they were used only with `Real` arguments (like all other mathematic functions) the ADL would find them.
// Also, the ADL does not work properly for g++ ver. older than 6 in case when Real is a fundamental type (e.g. `double`) inside `math::`
//    #if (__GNUC__ <= 6) …… #endif
using math::abs;
using math::fabs;
using math::max;
using math::min;

}

#undef YADE_WRAP_FUNC_1
#undef YADE_WRAP_FUNC_1_RENAME
#undef YADE_WRAP_FUNC_2
#undef YADE_WRAP_FUNC_2_TYPE2
#undef YADE_WRAP_FUNC_2_TYPE2_CAST
#undef YADE_WRAP_FUNC_2_TYPE2_STD_CAST
#undef YADE_WRAP_FUNC_3
#undef YADE_WRAP_FUNC_3_TYPE3
#undef YADE_WRAP_FUNC_1_COMPLEX
#undef YADE_WRAP_FUNC_1_COMPLEX_STD
#undef YADE_WRAP_FUNC_1_COMPLEX_TO_REAL
#undef YADE_WRAP_FUNC_1_COMPLEX_TO_REAL_STD

#endif

