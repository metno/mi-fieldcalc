/*
  mi-fieldcalc

  Copyright (C) 2019 met.no

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

#ifndef MI_FIELDCALC_OPENMP_TOOLS_H
#define MI_FIELDCALC_OPENMP_TOOLS_H 1

#define MIUTIL_STRING_HELPER0(x) #x
#define MIUTIL_STRING_HELPER1(x) MIUTIL_STRING_HELPER0(x)

#ifdef _OPENMP
#include <omp.h>

#define MIUTIL_OPENMP(options)                  \
  _Pragma(MIUTIL_STRING_HELPER1(omp options))

#define MIUTIL_OPENMP_PARALLEL(loopsize,options)                        \
  omp_set_dynamic(0);                                                   \
  omp_set_num_threads(miutil::compute_num_threads(loopsize));           \
  _Pragma(MIUTIL_STRING_HELPER1(omp parallel options))

#else // !_OPENMP

#define MIUTIL_OPENMP(options)                    /* nothing */
#define MIUTIL_OPENMP_PARALLEL(loopsize, options) /* nothing */

#endif // !_OPENMP

namespace miutil {

/*! A little utility function that dynamically computes the best number of openmp threads
  based on the number of iterations (nx*ny) for example. */
int compute_num_threads(long loop_size);

} // namespace miutil

#endif // MI_FIELDCALC_OPENMP_TOOLS_H
