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

IF (PYTHON_INSTALL_LIBDIR)
  MESSAGE (STATUS "Configured PYTHON_INSTALL_LIBDIR='${PYTHON_INSTALL_LIBDIR}'")
ELSE ()
  EXECUTE_PROCESS(OUTPUT_VARIABLE PYTHON_INSTALL_LIBDIR
    RESULT_VARIABLE PYTHON_INSTALL_EXIT_STATUS
    COMMAND ${Python_EXECUTABLE} -c "
import sys
if sys.version_info >= (3, 10):
  import sysconfig as s
  pl = s.get_path('platlib', 'posix_prefix', vars={'platbase':'${CMAKE_INSTALL_PREFIX}'})
else:
  from distutils import sysconfig as ds
  pl = ds.get_python_lib(plat_specific=True, prefix='${CMAKE_INSTALL_PREFIX}')
sys.stdout.write(pl)
")
  IF (PYTHON_INSTALL_EXIT_STATUS)
    MESSAGE(FATAL_ERROR "Could not detect PYTHON_INSTALL_LIBDIR")
  ELSE ()
    MESSAGE (STATUS "Detected PYTHON_INSTALL_LIBDIR='${PYTHON_INSTALL_LIBDIR}'")
  ENDIF()
ENDIF()

SET(py_SOURCES
  py_mi_fieldcalc.cc
)

SET(py_MODULE mi_fieldcalc)

PYBIND11_ADD_MODULE(${py_MODULE} ${py_SOURCES})
TARGET_LINK_LIBRARIES(${py_MODULE}
  PRIVATE
  mi-fieldcalc
)
# see https://pybind11.readthedocs.io/en/stable/faq.html#recursive-template-instantiation-exceeded-maximum-depth-of-256
# the culprits seem to be functions vesselIcingModStall and vesselIcingMincog which have 17 arguments
TARGET_COMPILE_OPTIONS(${py_MODULE}
  PRIVATE
  -ftemplate-depth=2048
)

INSTALL(TARGETS ${py_MODULE}
  LIBRARY DESTINATION "${PYTHON_INSTALL_LIBDIR}"
)

SET(py_TESTS
  test_mi_fieldcalc.py
)

IF(ENABLE_TESTS)
  FOREACH(T ${py_TESTS})
    ADD_TEST(NAME ${T}
      COMMAND ${Python_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/${T}"
      )
  ENDFOREACH()
  SET_TESTS_PROPERTIES(${py_TESTS}
    PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}"
    )
ENDIF()

# dummy target to tell IDEs that these files are part of the project
ADD_CUSTOM_TARGET(dummy_testfiles SOURCES ${py_TESTS})
