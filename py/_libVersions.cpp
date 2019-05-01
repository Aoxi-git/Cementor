#include <core/Omega.hpp>
#include <lib/pyutil/doc_opts.hpp>
#include <string>
#include <boost/lexical_cast.hpp>
#include <boost/python/tuple.hpp>

namespace py = boost::python;

// I used the list in https://yade-dem.org/doc/installation.html#prerequisites
// If we need a version of some library not listed in doc/sphinx/installation.rst, then it must also be added to that list!

// 1. compiler
// https://gcc.gnu.org/onlinedocs/cpp/Common-Predefined-Macros.html
#ifdef __GNUC__
	py::list gccVer() {
		py::list ret;
		ret.append( py::make_tuple( __GNUC__   , __GNUC_MINOR__ , __GNUC_PATCHLEVEL__ ));
		std::string prefix{};
#ifndef BOOST_GCC
// https://www.boost.org/doc/libs/1_44_0/libs/config/doc/html/boost_config/boost_macro_reference.html#boost_config.boost_macro_reference.boost_informational_macros
		// use boost macro to detect if it's really a gcc compiler.
		prefix += "not_gcc_";
#endif
#ifdef BOOST_INTEL
		prefix += "INTEL_";
#endif
#ifdef BOOST_CLANG
		prefix += "CLANG_";
#endif
		ret.append( prefix+boost::lexical_cast<std::string>(__GNUC__)+"."+boost::lexical_cast<std::string>(__GNUC_MINOR__)+"."+boost::lexical_cast<std::string>(__GNUC_PATCHLEVEL__ ));
		return ret;
	}
#else
	py::list gccVer() { return {}; }
#endif

// http://clang.llvm.org/docs/LanguageExtensions.html#builtin-macros
#ifdef __clang__
	py::list clangVer() {
		py::list ret;
		ret.append( py::make_tuple( __clang_major__ , __clang_minor__ , __clang_patchlevel__ ));
		ret.append( boost::lexical_cast<std::string>( __clang_version__ ));
		return ret;
	}
#else
	py::list clangVer() { return {}; }
#endif

// 2. Boost
// https://www.boost.org/doc/libs/1_44_0/libs/config/doc/html/boost_config/boost_macro_reference.html#boost_config.boost_macro_reference.boost_informational_macros
	py::list boostVer() {
		py::list ret;
		ret.append( py::make_tuple( BOOST_VERSION / 100000 , BOOST_VERSION / 100 % 1000 , BOOST_VERSION % 100 ));
// that's the default BOOST_VERSION
//		ret.append( boost::lexical_cast<std::string>( BOOST_VERSION ));
// But the string with dots is a bit more readable
		ret.append( boost::lexical_cast<std::string>(BOOST_VERSION / 100000)+"."+boost::lexical_cast<std::string>(BOOST_VERSION / 100 % 1000)+"."+boost::lexical_cast<std::string>(BOOST_VERSION % 100) );
		return ret;
	}

#ifdef YADE_MPI
	#include <mpi.h>
// https://www.open-mpi.org/software/ompi/versions/
	py::list mpiVer() {
		py::list ret;
		ret.append( py::make_tuple                  (OMPI_MAJOR_VERSION   ,                                   OMPI_MINOR_VERSION  ,                                    OMPI_RELEASE_VERSION ));
// I didn't find a way to obtain a version string for mpi, so I construct my own. But if there is some MPI version string, better put it here.
//		ret.append( boost::lexical_cast<std::string>(OMPI_VERSION) );
		ret.append( boost::lexical_cast<std::string>(OMPI_MAJOR_VERSION)+"."+boost::lexical_cast<std::string>(OMPI_MINOR_VERSION)+"."+boost::lexical_cast<std::string>(OMPI_RELEASE_VERSION) );
		return ret;
	}
#else
	py::list mpiVer() { return {}; }
#endif

#ifdef YADE_OPENMP
// lib/base/openmp-accu.hpp
// FIXME - add it to doc/sphinx/installation.rst
#endif

#ifdef YADE_VTK
	#include <vtkVersion.h>
	py::list vtkVer() {
		py::list ret;
		ret.append( py::make_tuple (VTK_MAJOR_VERSION , VTK_MINOR_VERSION , VTK_BUILD_VERSION ));
		ret.append( boost::lexical_cast<std::string>(VTK_VERSION) );
		return ret;
	}
#else
        py::list vtkVer() { return {}; }
#endif

#ifdef YADE_CGAL
	#include <CGAL/version_macros.h>
	py::list cgalVer() {
		py::list ret;
		ret.append( py::make_tuple( CGAL_VERSION_MAJOR , CGAL_VERSION_MINOR , CGAL_VERSION_PATCH ) );
		ret.append( boost::lexical_cast<std::string>(CGAL_VERSION_STR) );
		return ret;
	}
#else
	py::list cgalVer() { return {}; }
#endif

py::dict allVersionsCpp(){
	py::dict ret;
	ret["gcc"   ] = gccVer();
	ret["clang" ] = clangVer();
	ret["boost" ] = boostVer();
	ret["mpi"   ] = mpiVer();
	ret["vtk"   ] = vtkVer();
	ret["cgal"  ] = cgalVer();
	return ret;
}

BOOST_PYTHON_MODULE(_libVersions){
	YADE_SET_DOCSTRING_OPTS;
/*
 * Maybe exporting them all is not necessary. All data is aggregated inside allVersionsCpp() anyway.
 *
	py::def("gccVer"  , gccVer   , "Returns g++   compiler version with which yade was compiled.                                                       \n\n:return: list in format ``[ (major,minor,patch) , \"versionString\" ]``.");
	py::def("clangVer", clangVer , "Returns clang compiler version with which yade was compiled.                                                       \n\n:return: list in format ``[ (major,minor,patch) , \"versionString\" ]``.");
	py::def("boostVer", boostVer , "Returns boost  library version with which yade was compiled.                                                       \n\n:return: list in format ``[ (major,minor,patch) , \"versionString\" ]``.");
	py::def("mpiVer"  , mpiVer   , "Returns MPI    library version with which yade was compiled, as detected from ``#include                 <mpi.h>``.\n\n:return: list in format ``[ (major,minor,patch) , \"versionString\" ]``.");
	py::def("vtkVer"  , vtkVer   , "Returns VTK    library version with which yade was compiled, as detected from ``#include          <vtkVersion.h>``.\n\n:return: list in format ``[ (major,minor,patch) , \"versionString\" ]``.");
	py::def("cgalVer" , cgalVer  , "Returns CGAL   library version with which yade was compiled, as detected from ``#include <CGAL/version_macros.h>``.\n\n:return: list in format ``[ (major,minor,patch) , \"versionString\" ]``.");
*/

// We can use C++ string literal just like """ """ in python to write docstrings (see. https://en.cppreference.com/w/cpp/language/string_literal )
// The """ is a custom delimeter, we could use    R"RAW( instead, or any other delimeter. This delimeter decides what will be the termination delimeter.

	py::def("allVersionsCpp", allVersionsCpp, R"""(
This function returns library versions as discovered by C++ during compilation from all the ``#include`` headers. This can be useful in debugging to detect some library ``.so`` conflicts.

:return: dictionary in folowing format: ``{ "libName" : [ (major,minor,patch) , "versionString" ] }``

As an example the table below reflects with what libraries this documentation was compiled (here are only those detected by C++):

.. ipython::

	   In [1]: from yade._libVersions import *

	   In [1]: allVersionsCpp()

.. note:: Please add here C++ detection of other libraries when yade starts using them.

	)""");
}

