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

INCLUDE_DIRECTORIES(
  ${CMAKE_SOURCE_DIR}/src
  ${CMAKE_CURRENT_SOURCE_DIR}
)

ADD_DEFINITIONS(
  -DTEST_SRCDIR=\"${CMAKE_CURRENT_SOURCE_DIR}\"
  -DTEST_BUILDDIR=\"${CMAKE_CURRENT_BINARY_DIR}\"
)

LINK_LIBRARIES(
  ${GTEST_LIBRARY}
  ${GTEST_MAIN_LIBRARY}
)

########################################################################

IF(ENABLE_TESTS)
  ADD_EXECUTABLE(mi_fieldcalc_test
    FieldCalculationsTest.cc
    MetConstantsTest.cc

    gtestMain.cc
  )

  TARGET_LINK_LIBRARIES(mi_fieldcalc_test
    mi-fieldcalc
    Threads::Threads
  )

  ADD_TEST(NAME mi_fieldcalc_test
    COMMAND mi_fieldcalc_test --gtest_color=yes
  )
ENDIF()
