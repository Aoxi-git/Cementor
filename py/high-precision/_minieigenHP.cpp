/*************************************************************************
*  2009-2012 © Václav Šmilauer                                           *
*  2019        Janek Kozicki                                             *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

// When yade uses high-precision number as Real type the usual (old) 'import minieigen'
// has to be replaced with a minieigenHP library which uses exacly the same number of decimal places
// as yade is using everywhere else. Please note that such precision can be very arbitrary, because cpp_bin_float
// or mpfr take it as a compile-time argument. Hence such minieigenHP cannot be separately precompiled as a package.
// Though it could be precompiled for some special types such as boost::multiprecision::float128

#include <lib/base/Logging.hpp>
#include <lib/high-precision/Real.hpp>
#include <lib/high-precision/ToFromPythonConverter.hpp>
#include <lib/pyutil/doc_opts.hpp>
#include <iomanip>
#include <iostream>
#include <limits>
#include <py/high-precision/_ExposeStorageOrdering.hpp>
#include <sstream>

using namespace ::yade::MathEigenTypes;
#include <py/high-precision/minieigen/converters.hpp>
#include <py/high-precision/minieigen/expose.hpp>
#include <py/high-precision/minieigen/visitors.hpp>

CREATE_CPP_LOCAL_LOGGER("_minieigenHP.cpp")

BOOST_PYTHON_MODULE(_minieigenHP)
try {
	YADE_SET_DOCSTRING_OPTS;

	// arbitrary Real specific stuff: start
#if YADE_REAL_BIT > 64
	// these are needed only for high precision. The float and double are covered by default converters.
	ArbitraryComplex_from_python<Complex>();
	py::to_python_converter<Complex, ArbitraryComplex_to_python<Complex>>();

	ArbitraryReal_from_python<Real>();
	py::to_python_converter<Real, ArbitraryReal_to_python<Real>>();
#endif

	expose_storage_ordering();
	// arbitrary Real specific stuff: end

	py::scope().attr("__doc__") = "miniEigen is wrapper for a small part of the `Eigen <http://eigen.tuxfamily.org>`_ library. Refer to its documentation "
	                              "for details. All classes in this module support pickling.";

	expose_converters(); // in expose-converters.cpp
#ifndef EIGEN_DONT_ALIGN
	// when we use vectorization the Vector3r is AlignedVector3, so we need to register converter from plain old Vector3<Real> so that other functions can accept it as an argument
	custom_VectorAnyAny_from_sequence<Eigen::Matrix<Real, 3, 1>>();
	py::class_<Eigen::Matrix<Real, 3, 1>>(
	        "Vector3na",
	        "3-dimensional non-aligned float vector; same as :obj:`Vector3`, but with alignment (``Eigen::AlignedVector3``).\n\nSupported operations "
	        "(``f`` if a float/int, ``v`` is a Vector3): ``-v``, ``v+v``, ``v+=v``, ``v-v``, ``v-=v``, ``v*f``, ``f*v``, ``v*=f``, ``v/f``, ``v/=f``, "
	        "``v==v``, ``v!=v``, plus operations with ``Matrix3`` and ``Quaternion``.\n\nImplicit conversion from sequence (list, tuple, ...) of 3 "
	        "floats.\n\nStatic attributes: ``Zero``, ``Ones``, ``UnitX``, ``UnitY``, ``UnitZ``.",
	        py::init<>())
	        .def(VectorVisitor<Eigen::Matrix<Real, 3, 1>>());
#endif

	expose_vectors1();
	expose_vectors2();
	expose_matrices1(); // must come after vectors
	expose_matrices2(); // must come after vectors
	expose_complex1();
	expose_complex2();
	expose_quaternion();
	expose_boxes();

	// Requires -ldouble-conversion, but that function isn't used anywhere
	//	py::def("float2str",&doubleToShortest,(py::arg("f"),py::arg("pad")=0),"Return the shortest string representation of *f* which will is equal to *f* when converted back to float. This function is only useful in Python prior to 3.0; starting from that version, standard string conversion does just that.");

#ifdef EIGEN_DONT_ALIGN
	py::scope().attr("vectorize") = false;
#else
	py::scope().attr("vectorize") = true;
#endif
} catch (...) {
	LOG_FATAL("Importing this module caused an exception and this module is in an inconsistent state now.");
	PyErr_Print();
	PyErr_SetString(PyExc_SystemError, __FILE__);
	boost::python::handle_exception();
	throw;
}

