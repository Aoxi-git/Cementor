#ifdef _HIGH_PRECISION_SUPPORT
// compilation wall clock time: 0:24.23 → split into two files → 0:14.13
#include <lib/high-precision/Real.hpp>
#include <lib/high-precision/ToFromPythonConverter.hpp>
using namespace ::yade::MathEigenTypes;
// half of minieigen/expose-matrices.cpp
#include <py/high-precision/minieigen/visitors.hpp>
void expose_matrices2(){
	py::class_<MatrixXr>("MatrixX","XxX (dynamic-sized) float matrix. Constructed from list of rows (as VectorX).\n\nSupported operations (``m`` is a MatrixX, ``f`` if a float/int, ``v`` is a VectorX): ``-m``, ``m+m``, ``m+=m``, ``m-m``, ``m-=m``, ``m*f``, ``f*m``, ``m*=f``, ``m/f``, ``m/=f``, ``m*m``, ``m*=m``, ``m*v``, ``v*m``, ``m==m``, ``m!=m``.",py::init<>())
		.def(MatrixVisitor<MatrixXr>())
	;
}

#endif

