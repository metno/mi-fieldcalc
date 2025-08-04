/*
  mi-fieldcalc

  Copyright (C) 2019-2022 met.no

  Contact information:
  Norwegian Meteorological Institute
  Box 43 Blindern
  0313 OSLO
  NORWAY
  email: diana@met.no

  This file is part of mi-fieldcalc.

  mi-fieldcalc is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  mi-fieldcalc is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with mi-fieldcalc; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

#ifndef MI_FIELDCALC_VERSION_H
#define MI_FIELDCALC_VERSION_H

#define MI_FIELDCALC_VERSION_MAJOR 0
#define MI_FIELDCALC_VERSION_MINOR 1
#define MI_FIELDCALC_VERSION_PATCH 9

#define MI_FIELDCALC_VERSION_INT(major,minor,patch) \
    (1000000*major + 1000*minor + patch)
#define MI_FIELDCALC_VERSION_CURRENT_INT \
    MI_FIELDCALC_VERSION_INT(MI_FIELDCALC_VERSION_MAJOR, MI_FIELDCALC_VERSION_MINOR, MI_FIELDCALC_VERSION_PATCH)

#endif // MI_FIELDCALC_VERSION_H
