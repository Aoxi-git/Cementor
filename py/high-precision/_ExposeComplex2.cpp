/*************************************************************************
*  2012-2020 Václav Šmilauer                                             *
*  2020      Janek Kozicki                                               *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

// compilation wall clock time: 0:24.22 → split into two files → 0:16.20
#include <lib/high-precision/Real.hpp>
#include <lib/high-precision/ToFromPythonConverter.hpp>
using namespace ::yade::MathEigenTypes;
// half of minieigen/expose-complex.cpp file
#include <py/high-precision/minieigen/visitors.hpp>
template <int N> void expose_complex2(bool notDuplicate, const py::scope& topScope)
{
	if (notDuplicate) {
		py::class_<Vector6crHP<N>>("Vector6c", "/*TODO*/", py::init<>()).def(VectorVisitor<Vector6crHP<N>>());
	} else {
		py::scope().attr("Vector6c") = topScope.attr("Vector6c");
	}
}

// explicit instantination - tell compiler to produce a compiled version of expose_converters (it is faster when done in parallel in .cpp files)
YADE_HP_PYTHON_REGISTER(expose_complex2)
