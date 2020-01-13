// 2019 Janek Kozicki

/*
vtk changed their function name in vtk8, this is annoying. And we don't want to litter the code with
#ifdefs everywhere. Better to clean it up with vim comamnds:

:grep -E "InsertNextTupleValue" --include="*" . -R --exclude ChangeLog --exclude tags --exclude CMakeLists.txt
:%s/InsertNextTupleValue(\([^)]\+\))/INSERT_NEXT_TUPLE(\1)/gc

And use a macro in these places:
*/
#pragma once

#ifdef YADE_VTK
#include <vtkVersion.h>

// fix InsertNextTupleValue → InsertNextTuple name change

#if VTK_MAJOR_VERSION < 8
#define INSERT_NEXT_TUPLE(a) InsertNextTupleValue(a)
#define INSERT_NEXT_TYPED_TUPLE(a) InsertNextTupleValue(a)
#else
#define INSERT_NEXT_TUPLE(a) InsertNextTuple(a)
#define INSERT_NEXT_TYPED_TUPLE(a) InsertNextTypedTuple(a)
#endif

// it seems that VTK only accepts double, so let's convert it. Some precision is lost.
// But it is not calculations. It is only for recording of resuls.

#define INSERT_NEXT_POINT(val1, val2, val3) InsertNextPoint((static_cast<double>(val1)), (static_cast<double>(val2)), (static_cast<double>(val3)))

#define INSERT_NEXT_VALUE(val1) InsertNextValue((static_cast<double>(val1)))

#define ARRAY_3_DOUBLE(val1, val2, val3)                                                                                                                       \
	{                                                                                                                                                      \
		(static_cast<double>(val1)), (static_cast<double>(val2)), (static_cast<double>(val3))                                                          \
	}

#define ARRAY_9_DOUBLE(val1, val2, val3, val4, val5, val6, val7, val8, val9)                                                                                   \
	{                                                                                                                                                      \
		(static_cast<double>(val1)), (static_cast<double>(val2)), (static_cast<double>(val3)), (static_cast<double>(val4)),                            \
		        (static_cast<double>(val5)), (static_cast<double>(val6)), (static_cast<double>(val7)), (static_cast<double>(val8)),                    \
		        (static_cast<double>(val9)),                                                                                                           \
	}

// (and others in the future)
#endif


/*
# At first I wanted to do this inside cmake, but it turns out that function definitions
# are not supoprted: https://cmake.org/cmake/help/v3.0/prop_dir/COMPILE_DEFINITIONS.html#prop_dir:COMPILE_DEFINITIONS

IF(${VTK_MAJOR_VERSION} LESS 8)
  ADD_DEFINITIONS("-DINSERT_NEXT_TUPLE(a)=InsertNextTupleValue(a)")
ELSE()
  ADD_DEFINITIONS("-DINSERT_NEXT_TUPLE(a)=InsertNextTuple(a)")
ENDIF()
*/

