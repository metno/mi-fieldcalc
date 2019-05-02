# mi-fieldcalc
#
# Copyright (C) 2019 met.no
#
# Contact information:
# Norwegian Meteorological Institute
# Box 43 Blindern
# 0313 OSLO
# NORWAY
# email: diana@met.no
#
# This file is part of mi-fieldcalc.
#
# mi-fieldcalc is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# mi-fieldcalc is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with mi-fieldcalc; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

SET(mi_fieldcalc_SOURCES
  FieldCalculations.cc
  FieldDefined.cc
  MetConstants.cc
  openmp_tools.cc
)

METNO_HEADERS (mi_fieldcalc_HEADERS mi_fieldcalc_SOURCES ".cc" ".h")

LIST(APPEND mi_fieldcalc_HEADERS
  math_util.h
)

# these have no header file
LIST (APPEND mi_fieldcalc_SOURCES
  FieldCalculationsVesselIcing.cc
)

SET(mi_fieldcalc_SOURCES_OPT
  FieldCalculationsVesselIcing.cc
  FieldCalculations.cc
)

#IF (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
#  SET (OPT_FLAGS "${OPT_FLAGS} -Rpass-analysis=loop-vectorize")
#ENDIF ()

IF (CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  SET (OPT_FLAGS "${OPT_FLAGS} -mavx2")
ELSEIF (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
  SET (OPT_FLAGS "${OPT_FLAGS} -march=core-avx2")
ENDIF ()

SET_SOURCE_FILES_PROPERTIES(${mi_fieldcalc_SOURCES_OPT}
  PROPERTIES COMPILE_FLAGS "${OPT_FLAGS}"
)

########################################################################

ADD_LIBRARY(mi-fieldcalc SHARED
  ${mi_fieldcalc_SOURCES}
)

TARGET_INCLUDE_DIRECTORIES(mi-fieldcalc
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

SET_TARGET_PROPERTIES(mi-fieldcalc PROPERTIES
  VERSION "${MI_FIELDCALC_PVERSION_FULL}"
  SOVERSION "${MI_FIELDCALC_VERSION_MAJOR}"
)

INSTALL(TARGETS mi-fieldcalc
  EXPORT mi-fieldcalc
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

INSTALL(FILES
  ${mi_fieldcalc_HEADERS}
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/mi_fieldcalc"
)
