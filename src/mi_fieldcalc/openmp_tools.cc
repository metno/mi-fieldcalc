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

#include "openmp_tools.h"

#ifdef _OPENMP
#include <algorithm>
#endif

namespace miutil {

int compute_num_threads(long loop_size)
{
#ifdef _OPENMP
  /* We must look at these two wariables
     OMP_NUM_THREADS, OMP_THREAD_LIMIT
  */

  static int loc_omp_num_threads = -1;

  if (loc_omp_num_threads == -1) {
    if (getenv("OMP_NUM_THREADS") != 0) {
      loc_omp_num_threads = atoi(getenv("OMP_NUM_THREADS"));
    } else if (getenv("OMP_THREAD_LIMIT") != 0) {
      loc_omp_num_threads = atoi(getenv("OMP_THREAD_LIMIT"));
    } else {
      /* The default */
      loc_omp_num_threads = 8;
    }
  }

  if (loop_size < 1000)
    return std::min(1, loc_omp_num_threads);
  else if (loop_size <= 10000)
    return std::min(2, loc_omp_num_threads);
  else if (loop_size <= 100000)
    return std::min(4, loc_omp_num_threads);
  else
    return std::min(8, loc_omp_num_threads);
#else // !_OPENMP
  return 1;
#endif // !_OPENMP
}

} // namespace diutil
