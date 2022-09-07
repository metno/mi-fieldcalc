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

FIND_PACKAGE(GTest QUIET)
IF(GTEST_FOUND)
  IF (CMAKE_VERSION VERSION_LESS "3.20")
    SET(GTEST_LIBRARY GTest::GTest)
    SET(GTEST_MAIN_LIBRARY GTest::Main)
  ELSE()
    SET(GTEST_LIBRARY GTest::gtest)
    SET(GTEST_MAIN_LIBRARY GTest::gtest_main)
  ENDIF()
ELSE()
  MESSAGE("apparently no compiled GTest library, trying to build it")
  FIND_FILE(GTEST_DIR src/gtest-all.cc
    HINTS
    "${GTEST_ROOT}/src/gtest"
    /usr/src/googletest/googletest
    /usr/src/gmock/gtest
    /usr/src/gtest
    /usr/local/src/gtest
  )
  IF(NOT GTEST_DIR)
    MESSAGE(FATAL_ERROR "could not find gtest-all.cc")
  ENDIF()
  GET_FILENAME_COMPONENT(GTEST_DIR ${GTEST_DIR} PATH)
  GET_FILENAME_COMPONENT(GTEST_DIR ${GTEST_DIR} PATH)
  ADD_SUBDIRECTORY(${GTEST_DIR} ${CMAKE_CURRENT_BINARY_DIR}/gtest EXCLUDE_FROM_ALL)
  SET(GTEST_LIBRARY gtest)
  SET(GTEST_MAIN_LIBRARY gtest_main)
ENDIF()

