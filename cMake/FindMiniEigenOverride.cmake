##########################################################################
#  2019        Janek Kozicki                                             #
#                                                                        #
#  This program is free software; it is licensed under the terms of the  #
#  GNU General Public License v2 or later. See file LICENSE for details. #
##########################################################################
#
# This module sets following variable:
#
#  MINIEIGENOVERRIDE_FOUND - new version of python3-minieigen found


MESSAGE(STATUS "Checking file ${MINIEIGEN_PREFIX_PATH}/minieigen/visitors.hpp provided by newer version of package python3-minieigen")

IF((EXISTS "${MINIEIGEN_PREFIX_PATH}/minieigen/visitors.hpp") AND (EXISTS "${MINIEIGEN_PREFIX_PATH}/minieigen/common.hpp"))
	FILE(READ "${MINIEIGEN_PREFIX_PATH}/minieigen/common.hpp" COMMON_FILE)
	STRING(FIND "${COMMON_FILE}" "#ifndef MINIEIGEN_OVERRIDE" IS_OVERRIDE_PRESENT)
	IF(${IS_OVERRIDE_PRESENT} EQUAL -1)
		SET(MINIEIGENOVERRIDE_FOUND FALSE)
	ELSE()
		SET(MINIEIGENOVERRIDE_FOUND TRUE)
	ENDIF()
ELSE()
	SET(MINIEIGENOVERRIDE_FOUND FALSE)
ENDIF()

