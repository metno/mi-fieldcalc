/* -*- c++ -*-
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

#ifndef MI_FIELDCALC_MATH_UTIL_H
#define MI_FIELDCALC_MATH_UTIL_H 1

#include <algorithm>
#include <cmath>

namespace miutil {

template <typename T>
inline T square(T x)
{
  return x * x;
}

/**
 * Calculate the squared length of a 2-dim vector, x*x + y*y.
 */
template <typename T>
inline T absval2(T x, T y)
{
  return square(x) + square(y);
}

/**
 * Calculate the length of a 2-dim vector, sqrt(x*x + y*y).
 */
template <typename T>
inline T absval(T x, T y)
{
  return std::sqrt(absval2(x, y));
}

template <typename T>
inline void sort2(T& a, T& b)
{
  if (b < a)
    std::swap(a, b);
}

template <typename T1, typename T2>
void minimize(T1& a, const T2& b)
{
  if (b < a)
    a = b;
}

template <typename T1, typename T2>
void maximize(T1& a, const T2& b)
{
  if (b > a)
    a = b;
}

template <typename T1, typename T2>
void minimaximize(T1& mi, T1& ma, const T2& b)
{
  if (b > ma)
    ma = b;
  if (b < mi)
    mi = b;
}

template <typename T1>
bool value_between(const T1& v, const T1& lim0, const T1& lim1)
{
  if (lim0 <= lim1)
    return lim0 <= v and v <= lim1;
  else
    return lim1 <= v and v <= lim0;
}

template <typename T1>
T1 constrain_value(const T1& v, const T1& lim0, const T1& lim1)
{
  if (lim0 <= lim1) {
    if (v < lim0)
      return lim0;
    else if (lim1 < v)
      return lim1;
    else
      return v;
  } else {
    if (v < lim1)
      return lim1;
    else if (lim0 < v)
      return lim0;
    else
      return v;
  }
}

template <typename T>
T pow10(T t)
{
  return std::pow(10, t);
};

} // namespace miutil

#endif // MI_FIELDCALC_MATH_UTIL_H
