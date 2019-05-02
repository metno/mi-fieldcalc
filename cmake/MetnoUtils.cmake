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

FUNCTION(METNO_PVERSION_DEFINES pack header_file)
  FILE(STRINGS "${header_file}" version_definitions REGEX "^#define.*_VERSION_(MAJOR|MINOR|PATCH) +")
  FOREACH(v_def ${version_definitions})
    STRING(REGEX REPLACE "^#define.*_VERSION_(MAJOR|MINOR|PATCH) +([0-9]+) *$" \\1 v_type   "${v_def}")
    STRING(REGEX REPLACE "^#define.*_VERSION_(MAJOR|MINOR|PATCH) +([0-9]+) *$" \\2 v_number "${v_def}")
    SET(version_${v_type} "${v_number}")
  ENDFOREACH()

  SET(${pack}_VERSION_MAJOR "${version_MAJOR}" PARENT_SCOPE)
  SET(${pack}_VERSION_MINOR "${version_MINOR}" PARENT_SCOPE)
  SET(${pack}_VERSION_PATCH "${version_PATCH}" PARENT_SCOPE)

  SET(${pack}_PVERSION "${version_MAJOR}.${version_MINOR}" PARENT_SCOPE)
  SET(${pack}_PVERSION_FULL "${version_MAJOR}.${version_MINOR}.${version_PATCH}" PARENT_SCOPE)
ENDFUNCTION()

########################################################################

FUNCTION(METNO_HEADERS headers sources source_suffix header_suffix)
  FOREACH (src ${${sources}})
    STRING(REPLACE ${source_suffix} ${header_suffix} hdr ${src})
    LIST(APPEND lheaders ${hdr})
  ENDFOREACH ()
  SET(${headers} ${lheaders} PARENT_SCOPE)
ENDFUNCTION ()
