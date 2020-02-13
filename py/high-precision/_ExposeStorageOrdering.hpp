/*************************************************************************
*  2019 Janek Kozicki                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#include <boost/python.hpp>
#include <minieigen/common.hpp>

namespace py = ::boost::python;

template <typename T> int getEigenFlagTemplate() { return int(T::Flags); }

// extracts information about https://eigen.tuxfamily.org/dox/group__flags.html
py::dict getEigenFlags()
{
	py::dict ret {};
	//	ret["Index"]         = getEigenFlagTemplate<Index>(); // it's just typedef of `std::ptrdiff_t`, does not have options
	ret["Vector2i"]  = getEigenFlagTemplate<Vector2i>();
	ret["Vector2r"]  = getEigenFlagTemplate<Vector2r>();
	ret["Vector3i"]  = getEigenFlagTemplate<Vector3i>();
	ret["Vector3r"]  = getEigenFlagTemplate<Vector3r>();
	ret["Vector3ra"] = getEigenFlagTemplate<Vector3ra::Base>();
	ret["Vector4r"]  = getEigenFlagTemplate<Vector4r>();
	ret["Vector6i"]  = getEigenFlagTemplate<Vector6i>();
	ret["Vector6r"]  = getEigenFlagTemplate<Vector6r>();
	ret["Matrix3r"]  = getEigenFlagTemplate<Matrix3r>();
	ret["Matrix6r"]  = getEigenFlagTemplate<Matrix6r>();
	ret["MatrixXr"]  = getEigenFlagTemplate<MatrixXr>();
	ret["VectorXr"]  = getEigenFlagTemplate<VectorXr>();
	// /*
	ret["Vector2rDA"]    = getEigenFlagTemplate<Vector2rDA>();
	ret["Vector3rDA"]    = getEigenFlagTemplate<Vector3rDA>();
	ret["Vector4rDA"]    = getEigenFlagTemplate<Vector4rDA>();
	ret["Vector6rDA"]    = getEigenFlagTemplate<Vector6rDA>();
	ret["Matrix3rDA"]    = getEigenFlagTemplate<Matrix3rDA>();
	ret["Matrix6rDA"]    = getEigenFlagTemplate<Matrix6rDA>();
	ret["MatrixXrDA"]    = getEigenFlagTemplate<MatrixXrDA>();
	ret["VectorXrDA"]    = getEigenFlagTemplate<VectorXrDA>();
	ret["QuaternionrAL"] = getEigenFlagTemplate<QuaternionrAL>();
	// */
	ret["Quaternionr"]  = getEigenFlagTemplate<Quaternionr>();
	ret["AngleAxisr"]   = getEigenFlagTemplate<AngleAxisr::VectorType::Base>();
	ret["AlignedBox3r"] = getEigenFlagTemplate<AlignedBox3r::VectorType::Base>();
	ret["AlignedBox2r"] = getEigenFlagTemplate<AlignedBox2r::VectorType::Base>();
#ifdef _COMPLEX_SUPPORT
	ret["Vector2cr"] = getEigenFlagTemplate<Vector2cr>();
	ret["Vector3cr"] = getEigenFlagTemplate<Vector3cr>();
	ret["Vector6cr"] = getEigenFlagTemplate<Vector6cr>();
	ret["VectorXcr"] = getEigenFlagTemplate<VectorXcr>();
	ret["Matrix3cr"] = getEigenFlagTemplate<Matrix3cr>();
	ret["Matrix6cr"] = getEigenFlagTemplate<Matrix6cr>();
	ret["MatrixXcr"] = getEigenFlagTemplate<MatrixXcr>();
#endif
	return ret;
}

template <typename T> int getEigenStorageOrderTemplate() { return int(Eigen::internal::traits<T>::Options); }

// extracts information about https://eigen.tuxfamily.org/dox/group__TopicStorageOrders.html
py::dict getEigenStorageOrders()
{
	py::dict ret {};
	//	ret["Index"]         = getEigenStorageOrderTemplate<Index>(); // it's just typedef of `std::ptrdiff_t`, does not have options
	ret["Vector2i"]  = getEigenStorageOrderTemplate<Vector2i>();
	ret["Vector2r"]  = getEigenStorageOrderTemplate<Vector2r>();
	ret["Vector3i"]  = getEigenStorageOrderTemplate<Vector3i>();
	ret["Vector3r"]  = getEigenStorageOrderTemplate<Vector3r>();
	ret["Vector3ra"] = getEigenStorageOrderTemplate<Vector3ra>();
	ret["Vector4r"]  = getEigenStorageOrderTemplate<Vector4r>();
	ret["Vector6i"]  = getEigenStorageOrderTemplate<Vector6i>();
	ret["Vector6r"]  = getEigenStorageOrderTemplate<Vector6r>();
	ret["Matrix3r"]  = getEigenStorageOrderTemplate<Matrix3r>();
	ret["Matrix6r"]  = getEigenStorageOrderTemplate<Matrix6r>();
	ret["MatrixXr"]  = getEigenStorageOrderTemplate<MatrixXr>();
	ret["VectorXr"]  = getEigenStorageOrderTemplate<VectorXr>();
	// /*
	ret["Vector2rDA"]    = getEigenStorageOrderTemplate<Vector2rDA>();
	ret["Vector3rDA"]    = getEigenStorageOrderTemplate<Vector3rDA>();
	ret["Vector4rDA"]    = getEigenStorageOrderTemplate<Vector4rDA>();
	ret["Vector6rDA"]    = getEigenStorageOrderTemplate<Vector6rDA>();
	ret["Matrix3rDA"]    = getEigenStorageOrderTemplate<Matrix3rDA>();
	ret["Matrix6rDA"]    = getEigenStorageOrderTemplate<Matrix6rDA>();
	ret["MatrixXrDA"]    = getEigenStorageOrderTemplate<MatrixXrDA>();
	ret["VectorXrDA"]    = getEigenStorageOrderTemplate<VectorXrDA>();
	ret["QuaternionrAL"] = getEigenStorageOrderTemplate<QuaternionrAL::Coefficients>();
	// */
	ret["Quaternionr"]  = getEigenStorageOrderTemplate<Quaternionr::Coefficients>();
	ret["AngleAxisr"]   = getEigenStorageOrderTemplate<AngleAxisr::VectorType>();
	ret["AlignedBox3r"] = getEigenStorageOrderTemplate<AlignedBox3r::VectorType>();
	ret["AlignedBox2r"] = getEigenStorageOrderTemplate<AlignedBox2r::VectorType>();
#ifdef _COMPLEX_SUPPORT
	ret["Vector2cr"] = getEigenStorageOrderTemplate<Vector2cr>();
	ret["Vector3cr"] = getEigenStorageOrderTemplate<Vector3cr>();
	ret["Vector6cr"] = getEigenStorageOrderTemplate<Vector6cr>();
	ret["VectorXcr"] = getEigenStorageOrderTemplate<VectorXcr>();
	ret["Matrix3cr"] = getEigenStorageOrderTemplate<Matrix3cr>();
	ret["Matrix6cr"] = getEigenStorageOrderTemplate<Matrix6cr>();
	ret["MatrixXcr"] = getEigenStorageOrderTemplate<MatrixXcr>();
#endif
	return ret;
}

void expose_storage_ordering()
{
	py::def("getEigenStorageOrders", getEigenStorageOrders, R"""(
:return: A python dictionary listing options for all types, see: https://eigen.tuxfamily.org/dox/group__TopicStorageOrders.html
	)""");
	py::def("getEigenFlags", getEigenFlags, R"""(
:return: A python dictionary listing flags for all types, see: https://eigen.tuxfamily.org/dox/group__flags.html
	)""");
}

