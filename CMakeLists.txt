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

PROJECT(mi-fieldcalc LANGUAGES CXX)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.12)

SET(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
INCLUDE(MetnoUtils)
INCLUDE(CMakePackageConfigHelpers)
INCLUDE(GNUInstallDirs)

SET(CMAKE_CXX_STANDARD 11)

METNO_PVERSION_DEFINES(MI_FIELDCALC "src/mi_fieldcalc/mi_fieldcalc_version.h")

OPTION(ENABLE_OPENMP "Use OpenMP" OFF)
IF(ENABLE_OPENMP)
  FIND_PACKAGE(OpenMP)
  IF(OPENMP_FOUND)
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  ENDIF()
ENDIF()

SET(lib_name "mi-fieldcalc")

ENABLE_TESTING()

OPTION(ENABLE_TESTS "Enable tests" ON)
IF(ENABLE_TESTS)
  INCLUDE(MetnoGTest)
  SET(THREADS_PREFER_PTHREAD_FLAG TRUE)
  FIND_PACKAGE(Threads REQUIRED)
ENDIF()

ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(test)

OPTION(ENABLE_PYTHON "Enable python binding" ON)
IF(ENABLE_PYTHON)
  SET (PYBIND11_3RDPARTY "third_party/pybind11")
  IF ((NOT pybind11_DIR) AND (EXISTS "${CMAKE_SOURCE_DIR}/${PYBIND11_3RDPARTY}/CMakeLists.txt"))
    MESSAGE(STATUS "Including pybind11 project from ${PYBIND11_3RDPARTY}")
    ADD_SUBDIRECTORY(${PYBIND11_3RDPARTY})
  ELSE ()
    FIND_PACKAGE(pybind11 REQUIRED)
  ENDIF ()

  ADD_SUBDIRECTORY(python)
ENDIF()

SET(MI_FIELDCALC_CMAKE_DIR "${CMAKE_INSTALL_LIBDIR}/cmake/mi-fieldcalc")

INSTALL (
  EXPORT mi-fieldcalc
  DESTINATION "${MI_FIELDCALC_CMAKE_DIR}"
  FILE mi-fieldcalc-targets.cmake
)

WRITE_BASIC_PACKAGE_VERSION_FILE(
  "${CMAKE_CURRENT_BINARY_DIR}/mi-fieldcalc-config-version.cmake"
  VERSION
    "${MI_FIELDCALC_VERSION_MAJOR}.${MI_FIELDCALC_VERSION_MINOR}.${MI_FIELDCALC_VERSION_PATCH}"
  COMPATIBILITY
    SameMajorVersion
)

INSTALL(
  FILES
    mi-fieldcalc-config.cmake
    "${CMAKE_CURRENT_BINARY_DIR}/mi-fieldcalc-config-version.cmake"
  DESTINATION
    "${MI_FIELDCALC_CMAKE_DIR}"
)
