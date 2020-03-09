/*
 mi-fieldcalc - Meteorological Field Calculations

 Copyright (C) 2016-2019 met.no

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

#include "FieldCalculations.h"

#include "MetConstants.h"
#include "math_util.h"
#include "openmp_tools.h"

#include <memory>

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstring>

using namespace miutil::constants;

namespace {

template <typename T>
inline T pow2(T x)
{
  return x * x;
}

template <typename T>
inline T pow4(T x)
{
  return pow2(pow2(x));
}

void fillEdges(int nx, int ny, float* field)
{
  // Fill in edge values not computed or badly computed.
  // Note that in some (d/dx,d/dy) cases even correct computed edge values will be changed,
  // possibly not so good when current computation is used in following computations

  for (int j = 1, i1 = nx, i2 = nx * 2 - 1; j < ny - 1; j++, i1 += nx, i2 += nx) {
    field[i1] = field[i1 + 1];
    field[i2] = field[i2 - 1];
  }

  for (int i1 = 0, i2 = nx * ny - nx; i1 < nx; i1++, i2++) {
    field[i1] = field[i1 + nx];
    field[i2] = field[i2 - nx];
  }
}

bool fillUndef(int nx, int ny, float* fres, miutil::ValuesDefined& fDefined, float undef)
{
  const size_t fsize = nx * ny;
  fDefined = miutil::NONE_DEFINED;
  std::fill(fres, fres + fsize, undef);
  return true;
}

float divideUndef(float a, float b, float undef, size_t& n_undefined)
{
  if (b != 0) {
    return a / b;
  } else {
    n_undefined += 1;
    return undef;
  }
}

template <typename Func>
bool unaryFunctionField(const Func& func, int nx, int ny, const float* field, float* fres, miutil::ValuesDefined& fDefined, float undef)
{
  const size_t fsize = nx * ny;
  MIUTIL_OPENMP_PARALLEL(fsize, for)
  for (size_t i = 0; i < fsize; i++) {
    if (miutil::fieldcalc::is_defined(fDefined == miutil::ALL_DEFINED, field[i], undef))
      fres[i] = func(field[i]);
    else
      fres[i] = undef;
  }
  return true;
}

typedef float (*unary_function_t)(float f);

bool unaryFunctionField(unary_function_t func, int nx, int ny, const float* field, float* fres, miutil::ValuesDefined& fDefined, float undef)
{
  const size_t fsize = nx * ny;
  const bool inAllDefined = fDefined == miutil::ALL_DEFINED;
  MIUTIL_OPENMP_PARALLEL(fsize, for)
  for (size_t i = 0; i < fsize; i++) {
    if (miutil::fieldcalc::is_defined(inAllDefined, field[i], undef))
      fres[i] = func(field[i]);
    else
      fres[i] = undef;
  }
  return true;
}

typedef float (*binary_function_t)(float f1, float f2);

template <class Func>
bool binaryFunctionFieldField(const Func& func, int nx, int ny, const float* field1, const float* field2, float* fres, miutil::ValuesDefined& fDefined,
                              float undef)
{
  const size_t fsize = nx * ny;
  const bool inAllDefined = fDefined == miutil::ALL_DEFINED;
  MIUTIL_OPENMP_PARALLEL(fsize, for)
  for (size_t i = 0; i < fsize; i++) {
    if (miutil::fieldcalc::is_defined(inAllDefined, field1[i], field2[i], undef))
      fres[i] = func(field1[i], field2[i]);
    else
      fres[i] = undef;
  }
  return true;
}

template <typename Func>
bool unaryFunctionFieldUndef(const Func& func, int nx, int ny, const float* field, float* fres, miutil::ValuesDefined& fDefined, float undef)
{
  const size_t fsize = nx * ny;
  const bool inAllDefined = fDefined == miutil::ALL_DEFINED;
  size_t n_undefined = 0;
  MIUTIL_OPENMP_PARALLEL(fsize, for reduction(+:n_undefined))
  for (size_t i = 0; i < fsize; i++) {
    if (miutil::fieldcalc::is_defined(inAllDefined, field[i], undef)) {
      fres[i] = func(field[i], undef, n_undefined);
    } else {
      fres[i] = undef;
      n_undefined += 1;
    }
  }
  fDefined = miutil::checkDefined(n_undefined, fsize);
  return true;
}

template <class Func>
bool binaryFunctionFieldFieldUndef(const Func& func, int nx, int ny, const float* field1, const float* field2, float* fres, miutil::ValuesDefined& fDefined,
                                   float undef)
{
  const int fsize = nx * ny;
  const bool inAllDefined = fDefined == miutil::ALL_DEFINED;
  size_t n_undefined = 0;
  MIUTIL_OPENMP_PARALLEL(fsize, for reduction(+:n_undefined))
  for (int i = 0; i < fsize; i++) {
    if (miutil::fieldcalc::is_defined(inAllDefined, field1[i], field2[i], undef)) {
      fres[i] = func(field1[i], field2[i], undef, n_undefined);
    } else {
      fres[i] = undef;
      n_undefined += 1;
    }
  }
  fDefined = miutil::checkDefined(n_undefined, fsize);
  return true;
}

} // namespace

namespace miutil {
namespace fieldcalc {

inline float clamp_rh(float rh)
{
  if (rh < rhmin)
    return rhmin;
  else if (rh > rhmax)
    return rhmax;
  else
    return rh;
}

inline float t_thesat(float tk, float p, float pi, float undef, size_t& n_undefined)
{
  ewt_calculator ewt(tk - t0);
  if (not ewt.defined()) {
    n_undefined += 1;
    return undef;
  }
  const float qsat = eps * ewt.value() / p;
  return (cp * tk + xlh * qsat) / pi;
}

inline float th_thesat(float th, float p, float pi, float undef, size_t& n_undefined)
{
  ewt_calculator ewt(th * pi / cp - t0);
  if (not ewt.defined()) {
    n_undefined += 1;
    return undef;
  }
  const float qsat = eps * ewt.value() / p;
  return th + xlh * qsat / pi;
}

inline float tk_q_rh(float tk, float q, float p, float undef, size_t& n_undefined)
{
  ewt_calculator ewt(tk - t0);
  if (not ewt.defined()) {
    n_undefined += 1;
    return undef;
  }
  const float qsat = eps * ewt.value() / p;
  return 100. * q / qsat;
}

inline float tk_rh_q(float tk, float rh, float p, float undef, size_t& n_undefined)
{
  ewt_calculator ewt(tk - t0);
  if (not ewt.defined()) {
    n_undefined += 1;
    return undef;
  }
  const float qsat = eps * ewt.value() / p;
  return 0.01 * rh * qsat;
}

inline float tk_q_td(float tk, float q, float p, float tdconv, float undef, size_t& n_undefined)
{
  ewt_calculator ewt(tk - t0);
  if (not ewt.defined()) {
    n_undefined += 1;
    return undef;
  }

  const float et = ewt.value();
  const float qsat = eps * et / p;
  const float rh = fieldcalc::clamp_rh(q / qsat);
  const float etd = rh * et;
  return ewt.inverse(etd) + tdconv;
}

inline float tk_rh_td(float tk, float rh100, float tdconv, float undef, size_t& n_undefined)
{
  ewt_calculator ewt(tk - t0);
  if (not ewt.defined()) {
    n_undefined += 1;
    return undef;
  }

  const float et = ewt.value();
  const float rh = fieldcalc::clamp_rh(0.01 * rh100);
  const float etd = rh * et;
  return ewt.inverse(etd) + tdconv;
}

inline float tk_rh_the(float tk, float rh, float thconv, float undef, size_t& n_undefined)
{
  ewt_calculator ewt(tk - t0);
  if (not ewt.defined()) {
    n_undefined += 1;
    return undef;
  }

  return tk * thconv + ewt.value() * rh;
}

inline float tk_q_duct(float tk, float q, float p)
{
  return 77.6 * (p / tk) + 373000. * (q * p) / (eps * tk * tk);
}

inline float tk_rh_duct(float tk, float q, float p, float undef, size_t& n_undefined)
{
  ewt_calculator ewt(tk - t0);
  if (not ewt.defined()) {
    n_undefined += 1;
    return undef;
  }

  const float et = ewt.value();
  const float rh = clamp_rh(q * 0.01);
  return 77.6 * (p / tk) + 373000. * rh * et / (tk * tk);
}

float bad_hlevel(float a, float b)
{
  return (a < 0.0) || (b < 0.0) || (a == 0.0 && b == 0.0) || (b > 1.0);
}

inline float p_hlevel(float ps, float a, float b)
{
  return a + b * ps;
}

inline float pidcp_from_p(float p)
{
  return powf(p * p0inv, kappa);
}

inline float pi_from_p(float p)
{
  return cp * pidcp_from_p(p);
}

void copy_field(float* fout, const float* fin, size_t fsize)
{
  if (fout != fin)
    std::memcpy(fout, fin, sizeof(fin[0]) * fsize);
}

//---------------------------------------------------
// pressure level (PLEVEL) functions
//---------------------------------------------------

bool pleveltemp(int nx, int ny, const float* tinp, float p, const std::string& unit, int compute, float* tout, ValuesDefined& fDefined, float undef)
{
  if (p <= 0)
    return false;

  //  Pressure levels:
  //     compute=1 : potensiell temp. -> temp. (grader Celsius)
  //     compute=2 : potensiell temp. -> temp. (grader Kelvin)
  //     compute=3 : temp. (grader Kelvin)  -> potensiell temp. (Kelvin)
  //     compute=4 : temp. (grader Kelvin) -> saturated equivalent pot.temp.
  //     compute=5 : potensiell temp.      -> saturated equivalent pot.temp.

  if (compute < 3) {
    if (unit == "celsius")
      compute = 1;
    else if (unit == "kelvin")
      compute = 2;
  }

  const float pidcp = fieldcalc::pidcp_from_p(p), pi = pidcp * cp;

  switch (compute) {
  case 1: // TH -> T(Celsius)
    return unaryFunctionField([pidcp](float f) { return f * pidcp - t0; }, nx, ny, tinp, tout, fDefined, undef);
  case 2: // TH -> T(Kelvin)
    return unaryFunctionField([pidcp](float f) { return f * pidcp; }, nx, ny, tinp, tout, fDefined, undef);
  case 3: // T(Kelvin) -> TH
    return unaryFunctionField([pidcp](float f) { return f / pidcp; }, nx, ny, tinp, tout, fDefined, undef);
  case 4: { // T(Kelvin) -> THESAT
    auto func = [p, pi](float f, float undef, size_t& n_undefined) { return fieldcalc::t_thesat(f, p, pi, undef, n_undefined); };
    return unaryFunctionFieldUndef(func, nx, ny, tinp, tout, fDefined, undef);
  }
  case 5: { // TH -> THESAT
    auto func = [p, pi](float f, float undef, size_t& n_undefined) { return fieldcalc::th_thesat(f, p, pi, undef, n_undefined); };
    return unaryFunctionFieldUndef(func, nx, ny, tinp, tout, fDefined, undef);
  }
  default:
    return false;
  }
}

bool plevelthe(int nx, int ny, const float* t, const float* rh, float p, int compute, float* the, ValuesDefined& fDefined, float undef)
{
  //  Pressure levels:
  //    compute=1 : temp. (Kelvin)  og RH(%)     -> THE, ekvivalent pot.temp (Kelvin)
  //    compute=2 : pot. temp. (Kelvin) og RH(%) -> THE, ekvivalent pot.temp (Kelvin)
  //
  //  input:
  //    nx,ny      : field dimensions
  //    t[nx*ny]   : temp. (Kelvin)  or  pot.temp. (Kelvin)
  //    rh[nx*ny]  : rel. humidity (%)
  //    p          : pressure (hPa)
  //  output:
  //    the[nx*ny] : equivale pot. temp. (Kelvin)

  if (compute != 1 and compute != 2)
    return false;

  if (p <= 0.0)
    return false;

  const float pi = fieldcalc::pi_from_p(p);
  const float cvrh = 0.01 * (xlh / pi) * eps / p;
  const float tconv = (compute == 2) ? (pi / cp) : 1;
  const float thconv = (compute == 1) ? (cp / pi) : 1;

  auto func = [tconv, cvrh, thconv](float t, float rh, float undef, size_t& n_undefined) {
    return fieldcalc::tk_rh_the(t * tconv, rh * cvrh, thconv, undef, n_undefined);
  };
  return binaryFunctionFieldFieldUndef(func, nx, ny, t, rh, the, fDefined, undef);
}

bool plevelhum(int nx, int ny, const float* t, const float* huminp, float p,
               const std::string& unit, int compute, float* humout, ValuesDefined& fDefined, float undef)
{
  //  Pressure levels:
  //     compute=1  : temp. (Kelvin) og spes. fukt. -> rel. fuktighet (%)
  //     compute=2  : pot. temp. og spesifikk fukt. -> rel. fuktighet (%)
  //     compute=3  : temp. (Kelvin) og rel. fukt.  -> spesifikk fukt. (kg/kg)
  //     compute=4  : pot. temp. og  relativ fukt.  -> spesifikk fukt. (kg/kg)
  //     compute=5  : temp. (Kelvin) og rel. fukt.  -> duggpunkt, Td (Celsius)
  //     compute=6  : pot. temp. og  relativ fukt.  -> duggpunkt, Td (Celsius)
  //     compute=7  : temp. (Kelvin) og spes. fukt. -> duggpunkt, Td (Celsius)
  //     compute=8  : pot. temp. og spesifikk fukt. -> duggpunkt, Td (Celsius)
  //     compute=9  : temp. (Kelvin) og rel. fukt.  -> duggpunkt, Td (Kelvin)
  //     compute=10 : pot. temp. og  relativ fukt.  -> duggpunkt, Td (Kelvin)
  //     compute=11 : temp. (Kelvin) og spes. fukt. -> duggpunkt, Td (Kelvin)
  //     compute=12 : pot. temp. og spesifikk fukt. -> duggpunkt, Td (Kelvin)
  //
  //  Note: compute=5,9 are pressure independent (and not needed here)

  if (p <= 0 or compute <= 0 or compute >= 13)
    return false;

  if (compute > 8 && unit == "celsius")
    compute -= 4;
  else if (compute > 4 && compute <= 8 && unit == "kelvin")
    compute += 4;

  const int fsize = nx * ny;

  if (p == undef and (compute != 5 and compute != 6 and compute != 9 and compute != 10)) {
    // compute != 5/6/9/10 depends on p, so if p == undef, the result must also be undef
    return fillUndef(nx, ny, humout, fDefined, undef);
  }

  const float pi = fieldcalc::pi_from_p(p);

  const float tconv = (compute % 2 == 0) ? (pi / cp) : 1;
  const float tdconv = (compute >= 9) ? t0 : 0;

  const bool inAllDefined = fDefined == miutil::ALL_DEFINED;
  size_t n_undefined = 0;
  MIUTIL_OPENMP_PARALLEL(fsize, for reduction(+:n_undefined))
  for (int i = 0; i < fsize; i++) {
    if (fieldcalc::is_defined(inAllDefined, t[i], huminp[i], undef)) { // p checked before loop
      if (compute == 1 or compute == 2) {                                 // T(Kelvin),q -> RH(%)  or  TH,q -> RH(%)
        humout[i] = fieldcalc::tk_q_rh(t[i] * tconv, huminp[i], p, undef, n_undefined);
      } else if (compute == 3 or compute == 4) { // T(Kelvin),RH(%) -> q  or  TH,RH(%) -> q
        humout[i] = fieldcalc::tk_rh_q(t[i] * tconv, huminp[i], p, undef, n_undefined);
      } else if (compute == 5 or compute == 6 or compute == 9 or compute == 10) {
        // T(Kelvin),RH(%) -> Td(Celsius)     or  TH,RH(%) -> Td(Celsius)
        // or  T(Kelvin),RH(%) -> Td(Kelvin)  or  TH,RH(%) -> Td(Kelvin)
        humout[i] = fieldcalc::tk_rh_td(t[i] * tconv, huminp[i], tdconv, undef, n_undefined);
      } else if (compute == 7 or compute == 8 or compute == 11 or compute == 12) {
        // T(Kelvin),q -> Td(Celsius)     or  TH,q -> Td(Celsius)
        // or  T(Kelvin),q -> Td(Kelvin)  or  TH,q -> Td(Kelvin)
        humout[i] = fieldcalc::tk_q_td(t[i] * tconv, huminp[i], p, tdconv, undef, n_undefined);
      }
    } else {
      humout[i] = undef;
      n_undefined += 1;
    }
  }
  fDefined = miutil::checkDefined(n_undefined, fsize);
  return true;
}

bool pleveldz2tmean(int nx, int ny, const float* z1, const float* z2, float p1, float p2, int compute, float* tmean, ValuesDefined& fDefined, float undef)
{
  //  Pressure levels:
  //     compute=1 ; tykkelse -> middel temp. (grader Celsius)
  //     compute=2 ; tykkelse -> middel temp. (grader Kelvin)
  //     compute=3 ; tykkelse -> middel potensiell temp. (grader Kelvin)
  //  input: z1(nx,ny),z2(nx,ny) og p1,p2 (p1<p2 eller p1>p2)
  //  output: tmean(nx,ny) ........ tmean=z1 eller tmean=z2 er o.k.
  //
  //  antar middel-temp. gjelder i nivaa som er middel av exner-funk.

  if (p1 <= 0 || p2 <= 0 || p1 == p2)
    return false;

  const float pi1 = fieldcalc::pi_from_p(p1);
  const float pi2 = fieldcalc::pi_from_p(p2);

  float convert, tconvert;
  switch (compute) {
  case 1:
    convert = g * 0.5 * (pi1 + pi2) / ((pi2 - pi1) * cp);
    tconvert = -t0;
    break;
  case 2:
    convert = g * 0.5 * (pi1 + pi2) / ((pi2 - pi1) * cp);
    tconvert = 0.;
    break;
  case 3:
    convert = g / (pi2 - pi1);
    tconvert = 0.;
    break;
  default:
    return false;
  }

  auto func = [convert, tconvert](float z1, float z2) { return (z1 - z2) * convert + tconvert; };
  return binaryFunctionFieldField(func, nx, ny, z1, z2, tmean, fDefined, undef);
}

bool plevelqvector(int nx, int ny, const float* z, const float* t, const float* xmapr, const float* ymapr, const float* fcoriolis,
                   float p, int compute, float* qcomp, ValuesDefined& fDefined, float undef)
{
  //  Q-vector in pressure levels
  //
  //  compute=1 : input t is abs.temp. (K),       output qcomp (x-comp)
  //  compute=2 : input t is potential temp. (K), output qcomp (x-comp)
  //  compute=3 : input t is abs.temp. (K),       output qcomp (y-comp)
  //  compute=4 : input t is potential temp. (K), output qcomp (y-comp)
  //
  //  input:  z[nx*ny] : height (m) in pressure level
  //          t[nx*ny] : temperature or potential temp.
  //          xmapr[nx*ny],ymapr[nx*ny] : map ratios (xm/(hx*2),ym/(hy*2))
  //          fcoriolis[nx*ny] : coriolis parameter
  //          p: pressure (hPa)
  //  output: qx[nx*ny] : x component of Q-vector
  //          qy[nx*ny] : y component of Q-vector

  const int fsize = nx * ny;
  float tscale, pi;

  if (p <= 0.0)
    return false;

  if (nx < 3 || ny < 3)
    return false;

  switch (compute) {
  case 1:
    // t is abs.temp. (K)
    tscale = 1.0;
    break;
  case 2:
    // t is pot.temp. (K)
    pi = cp * powf(p / p0, r / cp);
    tscale = pi / cp;
    break;
  case 3:
    // t is abs.temp. (K)
    tscale = 1.0;
    break;
  case 4:
    // t is pot.temp. (K)
    pi = cp * powf(p / p0, r / cp);
    tscale = pi / cp;
    break;
  default:
    return false;
  }

  std::unique_ptr<float[]> ug(new float[fsize]);
  std::unique_ptr<float[]> vg(new float[fsize]);

  if (!plevelgwind_xcomp(nx, ny, z, xmapr, ymapr, fcoriolis, ug.get(), fDefined, undef) ||
      !plevelgwind_ycomp(nx, ny, z, xmapr, ymapr, fcoriolis, vg.get(), fDefined, undef))
  {
    return false;
  }

  const float c = -r / (p * 100.);

  // loop extended, reset bad computations at boundaries later
  size_t n_undefined = 0;
  MIUTIL_OPENMP_PARALLEL(fsize-2*nx, for reduction(+:n_undefined))
  for (int i = nx; i < fsize - nx; i++) {
    if (ug[i - nx] != undef && ug[i - 1] != undef && ug[i + 1] != undef && ug[i + nx] != undef && vg[i - nx] != undef && vg[i - 1] != undef &&
        vg[i + 1] != undef && vg[i + nx] != undef && t[i - nx] != undef && t[i - 1] != undef && t[i + 1] != undef && t[i + nx] != undef) {
      if (compute < 3) {
        const float dugdx = 0.5 * xmapr[i] * (ug[i + 1] - ug[i - 1]);
        const float dvgdx = 0.5 * xmapr[i] * (vg[i + 1] - vg[i - 1]);
        const float dtdx = 0.5 * xmapr[i] * tscale * (t[i + 1] - t[i - 1]);
        const float dtdy = 0.5 * ymapr[i] * tscale * (t[i + nx] - t[i - nx]);
        qcomp[i] = c * (dugdx * dtdx + dvgdx * dtdy);
      } else {
        const float dugdy = 0.5 * ymapr[i] * (ug[i + nx] - ug[i - nx]);
        const float dvgdy = 0.5 * ymapr[i] * (vg[i + nx] - vg[i - nx]);
        const float dtdx = 0.5 * xmapr[i] * tscale * (t[i + 1] - t[i - 1]);
        const float dtdy = 0.5 * ymapr[i] * tscale * (t[i + nx] - t[i - nx]);
        qcomp[i] = c * (dugdy * dtdx + dvgdy * dtdy);
      }
    } else {
      qcomp[i] = undef;
      n_undefined += 1;
    }
  }
  fDefined = miutil::checkDefined(n_undefined, fsize - 2 * nx);

  // fill in edge values not computed (or badly computed) above
  fillEdges(nx, ny, qcomp);
  return true;
}

bool plevelducting(int nx, int ny, const float* t, const float* h, float p, int compute, float* duct, ValuesDefined& fDefined, float undef)
{
  //  Pressure level:
  //
  //    compute=1 : temp. (Kelvin)      and q (kg/kg)  -> ducting
  //    compute=2 : pot. temp. (Kelvin) and q (kg/kg)  -> ducting
  //    compute=3 : temp. (Kelvin)      and RH(%)      -> ducting
  //    compute=4 : pot. temp. (Kelvin) and RH(%)      -> ducting
  //
  //  input:
  //    nx,ny      : field dimensions
  //    t[nx*ny]   : temp. (Kelvin)  or  pot.temp. (Kelvin)
  //    h[nx*ny]   : specific humidity (kg/kg) or rel.hum.(%)
  //    p          : pressure (hPa)
  //  output:
  //    duct[nx*ny]: ducting
  //...................................................................
  //       t in unit Kelvin, p in unit hPa
  //       duct=77.6*(p/t)+373000.*(q*p)/(eps*t*t)
  //       q*p/eps = rh*qsat*p/eps = rh*(eps*e(t)/p)*p/eps
  //               = rh*e(t) = (e(td)/e(t))*e(t) = e(td)
  //       => duct = 77.6*(p/t)+373000.*e(td)/(t*t)
  //          duct = 77.6*(p/t)+373000.*rh*e(t)/(t*t)
  //...................................................................

  if (p <= 0)
    return false;

  const float tconv = (compute % 2 == 0) ? fieldcalc::pidcp_from_p(p) : 1;

  if (compute == 1 or compute == 2) { // T(Kelvin),q -> ducting  or  TH,q -> ducting
    auto func = [tconv, p](float t, float h) { return fieldcalc::tk_q_duct(t * tconv, h, p); };
    return binaryFunctionFieldField(func, nx, ny, t, h, duct, fDefined, undef);
  } else if (compute == 3 or compute == 4) { // T(Kelvin),RH(%) -> ducting  or TH,RH(%) -> ducting
    auto func = [tconv, p](float t, float h, float undef, size_t& n_undefined) { return fieldcalc::tk_rh_duct(t * tconv, h, p, undef, n_undefined); };
    return binaryFunctionFieldFieldUndef(func, nx, ny, t, h, duct, fDefined, undef);
  } else {
    return false;
  }
}

bool plevelgwind_xcomp(int nx, int ny, const float* z, const float* /*xmapr*/, const float* ymapr, const float* fcoriolis, float* ug, ValuesDefined& fDefined,
                       float undef)
{
  //  Geostophic wind in pressure level
  //  (centered differences)
  //
  //  input:  z[nx*ny] : height (m) in pressure level
  //          xmapr[nx*ny],ymapr[nx*ny] : map ratios (xm/(hx*2),ym/(hy*2))
  //          fcoriolis[nx*ny] : coriolis parameter
  //  output: ug[nx*ny] : x component of geostrophic wind (m/s)

  if (nx < 3 || ny < 3)
    return false;

  const int fsize = nx * ny;

  // loop extended, reset bad computations at boundaries later
  const bool inAllDefined = fDefined == miutil::ALL_DEFINED;
  size_t n_undefined = 0;
  MIUTIL_OPENMP_PARALLEL(fsize-2*nx, for reduction(+:n_undefined))
  for (int i = nx; i < fsize - nx; i++) {
    // check for y component input, too
    if (fieldcalc::is_defined(inAllDefined, z[i - nx], z[i - 1], z[i + 1], z[i + nx], undef))
      ug[i] = -0.5 * ymapr[i] * (z[i + nx] - z[i - nx]) * g / fcoriolis[i];
    else
      ug[i] = undef;
    n_undefined += 1;
  }
  fDefined = miutil::checkDefined(n_undefined, fsize - 2 * nx);

  // fill in edge values not computed (or badly computed) above
  fillEdges(nx, ny, ug);

  return true;
}

bool plevelgwind_ycomp(int nx, int ny, const float* z, const float* xmapr, const float* /*ymapr*/, const float* fcoriolis, float* vg, ValuesDefined& fDefined,
                       float undef)
{
  //  Geostophic wind in pressure level
  //  (centered differences)
  //
  //  input:  z[nx*ny] : height (m) in pressure level
  //          xmapr[nx*ny],ymapr[nx*ny] : map ratios (xm/(hx*2),ym/(hy*2))
  //          fcoriolis[nx*ny] : coriolis parameter
  //  output: vg[nx*ny] : y component of geostrophic wind (m/s)

  const int fsize = nx * ny;

  // loop extended, reset bad computations at boundaries later
  const bool inAllDefined = fDefined == miutil::ALL_DEFINED;
  size_t n_undefined = 0;
  MIUTIL_OPENMP_PARALLEL(fsize-2*nx, for reduction(+:n_undefined))
  for (int i = nx; i < fsize - nx; i++) {
    // check for x component input, too
    if (fieldcalc::is_defined(inAllDefined, z[i - nx], z[i - 1], z[i + 1], z[i + nx], undef))
      vg[i] = 0.5 * xmapr[i] * (z[i + 1] - z[i - 1]) * g / fcoriolis[i];
    else {
      vg[i] = undef;
      n_undefined += 1;
    }
  }
  fDefined = miutil::checkDefined(n_undefined, fsize - 2 * nx);

  // fill in edge values not computed (or badly computed) above
  fillEdges(nx, ny, vg);

  return true;
}

bool plevelgvort(int nx, int ny, const float* z, const float* xmapr, const float* ymapr, const float* fcoriolis, float* gvort, ValuesDefined& fDefined,
                 float undef)
{
  //  Geostophic vorticity in pressure level (centered differences)
  //
  //  input:  z[nx*ny] : height (m)
  //          xmapr[nx*ny],ymapr[nx*ny] : map ratios (xm/(2*hx),ym/(2*hy) !!!)
  //          fcoriolis[nx*ny] : coriolis parameter
  //  output: gvort[nx*ny] : geostrophic relative vorticity

  const float g4 = g * 4.;

  int fsize = nx * ny;

  if (nx < 3 || ny < 3)
    return false;

  const bool inAllDefined = fDefined == miutil::ALL_DEFINED;
  size_t n_undefined = 0;
  MIUTIL_OPENMP_PARALLEL(fsize-2*nx, for reduction(+:n_undefined))
  for (int i = nx; i < fsize - nx; i++) {
    if (fieldcalc::is_defined(inAllDefined, z[i - nx], z[i - 1], z[i], z[i + 1], z[i + nx], undef))
      gvort[i] = (0.25 * xmapr[i] * xmapr[i] * (z[i - 1] - 2. * z[i] + z[i + 1]) + 0.25 * ymapr[i] * ymapr[i] * (z[i - nx] - 2. * z[i] + z[i + nx])) * g4 /
                 fcoriolis[i];
    else {
      gvort[i] = undef;
      n_undefined += 1;
    }
  }
  fDefined = miutil::checkDefined(n_undefined, fsize - 2 * nx);

  // fill in edge values not computed (or badly computed) above
  fillEdges(nx, ny, gvort);

  return true;
}

bool kIndex(int nx, int ny, const float* t500, const float* t700, const float* rh700, const float* t850, const float* rh850,
            float p500, float p700, float p850, int compute, float* kfield, ValuesDefined& fDefined, float undef)
{
  // K-index: (t+td)850 - (t-td)700 - (t)500
  //
  // note: levels are not fixed
  //
  //  Pressure levels:
  //    compute=1 : tPPP arrays contain temperature (Kelvin)
  //    compute=2 : tPPP arrays contain pot. temp. (Kelvin)
  //
  //  input:
  //    nx,ny        : field dimensions
  //     tPPP[nx*ny] : temp. (Kelvin)  or  pot.temp. (Kelvin)
  //    rhPPP[nx*ny] : rel. humidity (%)
  //    pPPP         : pressure (hPa)
  //  output:
  //    kfield[nx*ny] : K-index

  if (p500 <= 0.0 || p500 >= p700 || p700 >= p850)
    return false;

  float cvt500, cvt700, cvt850;
  switch (compute) {
  case 1: // Temp(Kelvin) input
    cvt500 = 1.;
    cvt700 = 1.;
    cvt850 = 1.;
    break;
  case 2: // Pot.temp. (TH) input
    cvt500 = fieldcalc::pidcp_from_p(p500);
    cvt700 = fieldcalc::pidcp_from_p(p700);
    cvt850 = fieldcalc::pidcp_from_p(p850);
    break;
  default:
    return false;
  }

  const int fsize = nx * ny;
  const bool inAllDefined = fDefined == miutil::ALL_DEFINED;
  size_t n_undefined = 0;
  MIUTIL_OPENMP_PARALLEL(fsize, for reduction(+:n_undefined))
  for (int i = 0; i < fsize; i++) {
    if (fieldcalc::is_defined(inAllDefined, t500[i], t700[i], rh700[i], t850[i], rh850[i], undef)) {
      // 850 hPa: rh,T -> Td ... rh*e(T) = e(Td) => Td = ?
      const float rh = fieldcalc::clamp_rh(0.01 * rh850[i]);
      const float tc850 = cvt850 * t850[i] - t0;
      const float tc700 = cvt700 * t700[i] - t0;
      const ewt_calculator ewt850(tc850), ewt700(tc700);
      if (not(ewt850.defined() and ewt700.defined())) {
        kfield[i] = undef;
        n_undefined += 1;
      } else {
        const float etd850 = ewt850.value() * rh;
        const float tdc850 = ewt850.inverse(etd850);
        // 700 hPa: rh,T -> Td ... rh*e(T) = e(Td) => Td = ?
        const float rh = fieldcalc::clamp_rh(0.01 * rh700[i]);
        const float etd700 = ewt700.value() * rh;
        const float tdc700 = ewt700.inverse(etd700);
        const float tc500 = cvt500 * t500[i] - t0;
        kfield[i] = (tc850 + tdc850) - (tc700 - tdc700) - tc500;
      }
    } else {
      kfield[i] = undef;
      n_undefined += 1;
    }
  }
  fDefined = miutil::checkDefined(n_undefined, fsize);
  return true;
}

bool ductingIndex(int nx, int ny, const float* t850, const float* rh850, float p850, int compute, float* duct, ValuesDefined& fDefined, float undef)
{
  //  Pressure levels:
  //     compute=1 ; temp. (Kelvin) og rel. fukt. (%)  -> ducting
  //     compute=2 ; pot. temp. og  relativ fukt. (%)  -> ducting
  //
  // ducting = nw(T) - nw(Td),   nw(t) = b*ew(t)/(t*t)
  //
  // NB! nivaaet er ikke fast, men normalt 850 hPa.

  const float bduct = 3.8e+5;

  if (p850 <= 0.0)
    return false;

  float tconvert;
  switch (compute) {
  case 1: // T(Kelvin),RH(%) input
    tconvert = 1.;
    break;
  case 2: // TH,RH(%) input
    tconvert = fieldcalc::pidcp_from_p(p850);
    break;
  default:
    return false;
  }

  const int fsize = nx * ny;
  const bool inAllDefined = fDefined == miutil::ALL_DEFINED;

  size_t n_undefined = 0;
  MIUTIL_OPENMP_PARALLEL(fsize, for reduction(+:n_undefined))
  for (int i = 0; i < fsize; i++) {
    if (fieldcalc::is_defined(inAllDefined, t850[i], rh850[i], undef)) {
      // 850 hPa: rh,T -> Td ... rh*e(T) = e(Td)
      const float rh = fieldcalc::clamp_rh(0.01 * rh850[i]);
      const float tk = t850[i] * tconvert;
      const ewt_calculator ewt(tk - t0);
      if (not ewt.defined()) {
        duct[i] = undef;
        n_undefined += 1;
      } else {
        const float et = ewt.value();
        const float etd = et * rh;
        const float tdk = ewt.inverse(etd) + t0;
        duct[i] = bduct * (et / (tk * tk) - etd / (tdk * tdk));
      }
    } else {
      duct[i] = undef;
      n_undefined += 1;
    }
  }
  fDefined = miutil::checkDefined(n_undefined, fsize);
  return true;
}

bool showalterIndex(int nx, int ny, const float* t500, const float* t850, const float* rh850, float p500, float p850, int compute, float* sfield, ValuesDefined& fDefined, float undef)
{
  // Showalter index:  t500 - tx500,
  //                   tx500 er "t850" hevet toerr-adiabatisk
  //                   til metning og videre langs
  //                   fuktigadiabat til 500 mb.
  //
  // metode: Benytter vanlig (modell-) metode for beregning av
  //         fuktigadiabat, dvs. hever foerst langs toerr-adiabat og
  //         justerer fuktighet og varme i oeverste nivaa ved noen
  //         iterasjoner. Fuktigheten som justeres er her ikke
  //         metnings-verdien nederst. (beregner ikke p,t for
  //         nivaa hvor metning inntreffer.)
  //
  //    compute=1 : tPPP arrays contain temperature (Kelvin)
  //    compute=2 : tPPP arrays contain pot. temp. (Kelvin)
  //
  //  input:
  //    nx,ny        : field dimensions
  //     tPPP[nx*ny] : temp. (Kelvin)  or  pot.temp. (Kelvin)
  //    rhPPP[nx*ny] : rel. humidity (%)
  //    pPPP         : pressure (hPa)
  //  output:
  //    sfield[nx*ny] : Showalter index
  //
  // note: levels are not fixed

  // antall iterasjoner for beregning av fuktigadiabat
  // (relativt mange p.g.a. 850 - 500 hPa uten "mellom-nivaa")

  if (p500 <= 0.0 || p500 >= p850)
    return false;

  const float pi500 = fieldcalc::pi_from_p(p500);
  const float pi850 = fieldcalc::pi_from_p(p850);

  float cvt500, cvt850, dryadiabat;
  switch (compute) {
  case 1:
    cvt500 = 1.;
    cvt850 = 1.;
    dryadiabat = cp * (cp / pi850) * (pi500 / cp);
    break;
  case 2:
    cvt500 = pi500 / cp;
    cvt850 = pi850 / cp;
    dryadiabat = cp * (pi500 / cp);
    break;
  default:
    return false;
  }

  const int niter = 7;
  const int fsize = nx * ny;

  const bool inAllDefined = fDefined == miutil::ALL_DEFINED;
  size_t n_undefined = 0;
  MIUTIL_OPENMP_PARALLEL(fsize, for reduction(+:n_undefined))
  for (int i = 0; i < fsize; i++) {
    if (fieldcalc::is_defined(inAllDefined, t500[i], t850[i], rh850[i], undef)) {
      const float tk500 = cvt500 * t500[i];
      const float tk850 = cvt850 * t850[i];
      const float rh = fieldcalc::clamp_rh(0.01 * rh850[i]);
      const ewt_calculator ewt(tk850 - t0);
      if (not ewt.defined()) {
        sfield[i] = undef;
        n_undefined += 1;
      } else {
        const float etd = ewt.value() * rh;

        // fuktigadiabat  (tcl: cp * grader kelvin)
        // lift preliminary along dry-adiabat (pot.temp. constant)
        float tcl = dryadiabat * t850[i];
        float qcl = eps * etd / p850;

        // adjust humidity and moisture in 'niter' iterations
        for (int n = 0; n < niter; n++) {
          const ewt_calculator ewt2(tcl / cp - t0);
          if (not ewt2.defined())
            break;
          const float esat = ewt2.value();
          const float qsat = eps * esat / p500;
          float dq = qcl - qsat;
          const float a1 = cplr * qcl / tcl;
          const float a2 = exl / tcl;
          dq = dq / (1. + a1 * a2);
          qcl = qcl - dq;
          tcl = tcl + dq * xlh;
        }

        const float tx500 = tcl / cp;
        sfield[i] = tk500 - tx500;
      }
    } else {
      n_undefined += 1;
    }
  }
  fDefined = miutil::checkDefined(n_undefined, fsize);
  return true;
}

bool boydenIndex(int nx, int ny, const float* t700, const float* z700, const float* z1000, float p700, float p1000, int compute, float* bfield, ValuesDefined& fDefined, float undef)
{
  // Boyden index:  (Z700-Z1000)/10 - Tc700 - 200.
  //
  //    compute=1 : tPPP arrays contain temperature (Kelvin)
  //    compute=2 : tPPP arrays contain pot. temp. (Kelvin)
  //
  //  input:
  //    nx,ny        : field dimensions
  //    tPPP[nx*ny]  : temp. (Kelvin)  or  pot.temp. (Kelvin)
  //    zPPP[nx*ny]  : height (m)
  //    pPPP         : pressure (hPa)
  //  output:
  //    bfield[nx*ny] : Boyden index
  //
  // note: levels are not fixed

  if (compute <= 0 or compute >= 3)
    return false;

  if (p700 <= 0.0 || p700 >= p1000)
    return false;

  const float pi700 = cp * powf(p700 / p0, r / cp);
  const float tconv = (compute == 2) ? pi700 / cp : 1;
  const int fsize = nx * ny;

  const bool inAllDefined = fDefined == miutil::ALL_DEFINED;
  size_t n_undefined = 0;
  MIUTIL_OPENMP_PARALLEL(fsize, for reduction(+:n_undefined))
  for (int i = 0; i < fsize; i++) {
    if (fieldcalc::is_defined(inAllDefined, t700[i], z700[i], z1000[i], undef)) {
      const float tc700 = t700[i] * tconv - t0;
      bfield[i] = (z700[i] - z1000[i]) / 10. - tc700 - 200.;
    } else {
      bfield[i] = undef;
      n_undefined += 1;
    }
  }
  fDefined = miutil::checkDefined(n_undefined, fsize);
  return true;
}

bool sweatIndex(int nx, int ny, const float* t850, const float* t500, const float* td850, const float* td500, const float* u850, const float* v850,
                const float* u500, const float* v500, float* sindex, ValuesDefined& fDefined, float undef)
{
  // Severe Weather Threat Index
  // Sweat index:12Td850 + 20(TTI - 49) + 2ff850 + ff500 + 125(sin(d500 -d850)) +0.2
  // TTI = t850+td850-2*t500
  const int fsize = nx * ny;

  const bool inAllDefined = fDefined == miutil::ALL_DEFINED;
  size_t n_undefined = 0;
  MIUTIL_OPENMP_PARALLEL(fsize, for reduction(+:n_undefined))
  for (int i = 0; i < fsize; i++) {
    if (fieldcalc::is_defined(inAllDefined, t850[i], t500[i], td850[i], td500[i], u850[i], v850[i], u500[i], v500[i], undef)) {
      float ff850 = miutil::absval(u850[i], v850[i]);
      float ff500 = miutil::absval(u500[i], v500[i]);
      float sind500_d850 = (u500[i] * v850[i] - v500[i] * u850[i]) / (ff850 * ff500);
      sindex[i] = 32 * td850[i] + 20 * t850[i] - 40 * t500[i] - 20 * 49 + 2 * miutil::ms2knots(ff850) + miutil::ms2knots(ff500) + 125 * (sind500_d850 + 0.2);
    } else {
      sindex[i] = undef;
      n_undefined += 1;
    }
  }
  fDefined = miutil::checkDefined(n_undefined, fsize);
  return true;
}

//---------------------------------------------------
// hybrid model level (HLEVEL) functions
//---------------------------------------------------

bool hleveltemp(int nx, int ny, const float* tinp, const float* ps, float alevel, float blevel,
                const std::string& unit, int compute, float* tout, ValuesDefined& fDefined, float undef)
{
  //  Modell-flater, sigma/eta(hybrid):
  //    sigma: alevel=ptop*(1.-sigma)   blevel=sigma  (p=pt+(ps[]-pt)*sigma)
  //    eta:   alevel,blevel
  //    p = alevel + blevel * ps[]
  //
  //     compute=1 : potensiell temp. -> temp. (grader Celsius)
  //     compute=2 : potensiell temp. -> temp. (grader Kelvin)
  //     compute=3 : temp. (grader Kelvin)  -> potensiell temp. (Kelvin)
  //     compute=4 : temp. (grader Kelvin) -> saturated equivalent pot.temp.
  //     compute=5 : potensiell temp.      -> saturated equivalent pot.temp.

  if (compute < 3) {
    if (unit == "celsius")
      compute = 1;
    else if (unit == "kelvin")
      compute = 2;
  }

  const int fsize = nx * ny;
  const bool inAllDefined = fDefined == miutil::ALL_DEFINED;

  if (fieldcalc::bad_hlevel(alevel, blevel)) {
    return false;
  }

  size_t n_undefined = 0;
  MIUTIL_OPENMP_PARALLEL(fsize, for reduction(+:n_undefined))
  for (int i = 0; i < fsize; i++) {
    if (fieldcalc::is_defined(inAllDefined, tinp[i], ps[i], undef)) {
      const float p = fieldcalc::p_hlevel(ps[i], alevel, blevel);
      const float pidcp = fieldcalc::pidcp_from_p(p);
      if (compute == 1) { // TH -> T(Celsius)
        tout[i] = tinp[i] * pidcp - t0;
      } else if (compute == 2) { // TH -> T(Kelvin)
        tout[i] = tinp[i] * pidcp;
      } else if (compute == 3) { // T(Kelvin) -> TH
        tout[i] = tinp[i] / pidcp;
      } else if (compute == 4) { // T(Kelvin) -> THESAT
        tout[i] = fieldcalc::t_thesat(tinp[i], p, pidcp * cp, undef, n_undefined);
      } else if (compute == 5) { // TH -> THESAT
        tout[i] = fieldcalc::th_thesat(tinp[i], p, pidcp * cp, undef, n_undefined);
      }
    } else {
      tout[i] = undef;
      n_undefined += 1;
    }
  }
  fDefined = miutil::checkDefined(n_undefined, fsize);
  return true;
}

bool hlevelthe(int nx, int ny, const float* t, const float* q, const float* ps, float alevel, float blevel, int compute, float* the, ValuesDefined& fDefined,
               float undef)
{
  //  Modell-flater, sigma/eta(hybrid):
  //    sigma: alevel=ptop*(1.-sigma)   blevel=sigma  (p=pt+(ps[]-pt)*sigma)
  //    eta:   alevel,blevel
  //       p = alevel + blevel * ps[]
  //
  //    compute=1 ; temp. (Kelvin)  og q     -> THE, ekvivalent pot.temp (Kelvin)
  //    compute=2 ; pot. temp. (Kelvin) og q -> THE, ekvivalent pot.temp (Kelvin)
  //
  //  input:
  //    nx,ny      : field dimensions
  //    t[nx*ny]   : temp. (Kelvin or Celsius)  or  pot.temp. (Kelvin)
  //    q[nx*ny]   : specific humidity (kg/kg)
  //    ps[nx*ny]  : surface pressure (hPa)
  //  output:
  //    the[nx*ny] : equivale pot. temp. (Kelvin)

  const int fsize = nx * ny;

  if (fieldcalc::bad_hlevel(alevel, blevel)) {
    return false;
  }

  const bool inAllDefined = fDefined == miutil::ALL_DEFINED;
  size_t n_undefined = 0;
  MIUTIL_OPENMP_PARALLEL(fsize, for reduction(+:n_undefined))
  for (int i = 0; i < fsize; i++) {
    if (fieldcalc::is_defined(inAllDefined, t[i], q[i], ps[i], undef)) {
      const float p = fieldcalc::p_hlevel(ps[i], alevel, blevel);
      const float pi = fieldcalc::pi_from_p(p);
      if (compute == 1) // T(Kelvin), q -> THE
        the[i] = (t[i] * cp + q[i] * xlh) / pi;
      else if (compute == 2) // TH, q -> THE
        the[i] = t[i] + q[i] * xlh / pi;
    } else {
      the[i] = undef;
      n_undefined += 1;
    }
  }
  fDefined = miutil::checkDefined(n_undefined, fsize);
  return true;
}

bool hlevelhum(int nx, int ny, const float* t, const float* huminp, const float* ps, float alevel, float blevel, const std::string& unit, int compute, float* humout,
               ValuesDefined& fDefined, float undef)
{
  //  Modell-flater, sigma/eta(hybrid):
  //    sigma: alevel=ptop*(1.-sigma)   blevel=sigma  (p=pt+(ps[]-pt)*sigma)
  //    eta:   alevel,blevel
  //       p = alevel + blevel * ps[]
  //
  //     compute=1  : temp. (Kelvin) og spes. fukt. -> rel. fuktighet (%)
  //     compute=2  : pot. temp. og spesifikk fukt. -> rel. fuktighet (%)
  //     compute=3  : temp. (Kelvin) og rel. fukt.  -> spesifikk fukt. (kg/kg)
  //     compute=4  : pot. temp. og  relativ fukt.  -> spesifikk fukt. (kg/kg)
  //     compute=5  : temp. (Kelvin) og spes. fukt. -> duggpunkt, Td (Celsius)
  //     compute=6  : pot. temp. og spesifikk fukt. -> duggpunkt, Td (Celsius)
  //     compute=7  : temp. (Kelvin) og rel. fukt.  -> duggpunkt, Td (Celsius)
  //     compute=8  : pot. temp. og  relativ fukt.  -> duggpunkt, Td (Celsius)
  //     compute=9  : temp. (Kelvin) og spes. fukt. -> duggpunkt, Td (Kelvin)
  //     compute=10 : pot. temp. og spesifikk fukt. -> duggpunkt, Td (Kelvin)
  //     compute=11 : temp. (Kelvin) og rel. fukt.  -> duggpunkt, Td (Kelvin)
  //     compute=12 : pot. temp. og  relativ fukt.  -> duggpunkt, Td (Kelvin)
  //
  //  Note: compute=7,11 are pressure independant (and not needed here)

  if (compute <= 0 or compute >= 13)
    return false;
  if (fieldcalc::bad_hlevel(alevel, blevel)) {
    return false;
  }

  if (compute > 8 && unit == "celsius")
    compute -= 4;
  else if (compute > 4 && compute <= 8 && unit == "kelvin")
    compute += 4;

  const int fsize = nx * ny;
  const bool inAllDefined = fDefined == miutil::ALL_DEFINED;
  const float tdconv = (compute >= 9) ? t0 : 0;
  const bool need_p = not(compute == 7 or compute == 11);

  size_t n_undefined = 0;
  MIUTIL_OPENMP_PARALLEL(fsize, for reduction(+:n_undefined))
  for (int i = 0; i < fsize; i++) {
    if (fieldcalc::is_defined(inAllDefined, t[i], huminp[i], undef) and ((not need_p) or inAllDefined or ps[i] != undef)) {
      const float p = need_p ? fieldcalc::p_hlevel(ps[i], alevel, blevel) : 0;
      if (compute == 1) { // T(Kelvin),q -> RH(%)
        humout[i] = fieldcalc::tk_q_rh(t[i], huminp[i], p, undef, n_undefined);
      } else if (compute == 2) { // TH,q -> RH(%)
        const float tk = t[i] * fieldcalc::pidcp_from_p(p);
        humout[i] = fieldcalc::tk_q_rh(tk, huminp[i], p, undef, n_undefined);
      } else if (compute == 3) { // T(Kelvin),RH(%) -> q
        humout[i] = fieldcalc::tk_rh_q(t[i], huminp[i], p, undef, n_undefined);
      } else if (compute == 4) { // TH,RH(%) -> q
        const float tk = t[i] * fieldcalc::pidcp_from_p(p);
        humout[i] = fieldcalc::tk_rh_q(tk, huminp[i], p, undef, n_undefined);
      } else if (compute == 5 or compute == 9) { // T(Kelvin),q -> Td(Celsius) or oT(Kelvin),q -> Td(Kelvin)
        humout[i] = fieldcalc::tk_q_td(t[i], huminp[i], p, tdconv, undef, n_undefined);
      } else if (compute == 6 or compute == 10) { // TH,q -> Td(Celsius)  or  TH,q -> Td(Kelvin)
        const float tk = t[i] * fieldcalc::pidcp_from_p(p);
        humout[i] = fieldcalc::tk_q_td(tk, huminp[i], p, tdconv, undef, n_undefined);
      } else if (compute == 7 or compute == 11) { // T(Kelvin),RH(%) -> Td(Celsius)  or  T(Kelvin),RH(%) -> Td(Kelvin)
        humout[i] = fieldcalc::tk_rh_td(t[i], huminp[i], tdconv, undef, n_undefined);
      } else if (compute == 8 or compute == 12) { // TH,RH(%) -> Td(Celsius)  or  TH,RH(%) -> Td(Kelvin)
        const float tk = t[i] * fieldcalc::pidcp_from_p(p);
        humout[i] = fieldcalc::tk_rh_td(tk, huminp[i], tdconv, undef, n_undefined);
      }
    } else {
      humout[i] = undef;
      n_undefined += 1;
    }
  }
  fDefined = miutil::checkDefined(n_undefined, fsize);
  return true;
}

bool hlevelducting(int nx, int ny, const float* t, const float* h, const float* ps, float alevel, float blevel, int compute, float* duct,
                   ValuesDefined& fDefined, float undef)
{
  //  Modell-flater, sigma/eta(hybrid):
  //    sigma: alevel=ptop*(1.-sigma)   blevel=sigma  (p=pt+(ps[]-pt)*sigma)
  //    eta:   alevel,blevel
  //       p = alevel + blevel * ps[]
  //
  //    compute=1 ; temp. (Kelvin)      and q (kg/kg)  -> ducting
  //    compute=2 ; pot. temp. (Kelvin) and q (kg/kg)  -> ducting
  //    compute=3 ; temp. (Kelvin)      and RH(%)      -> ducting
  //    compute=4 ; pot. temp. (Kelvin) and RH(%)      -> ducting
  //
  //  input:
  //    nx,ny      : field dimensions
  //    t[nx*ny]   : temp. (Kelvin)  or  pot.temp. (Kelvin)
  //    h[nx*ny]   : specific humidity (kg/kg) or rel.humidity (%)
  //    ps[nx*ny]  : surface pressure (hPa)
  //  output:
  //    duct[nx*ny]: ducting
  //...................................................................
  //       t in unit Kelvin, p in unit hPa ???
  //       duct=77.6*(p/t)+373000.*(q*p)/(eps*t*t)
  //       q*p/eps = rh*qsat*p/eps = rh*(eps*e(t)/p)*p/eps
  //               = rh*e(t) = (e(td)/e(t))*e(t) = e(td)
  //       => duct = 77.6*(p/t)+373000.*e(td)/(t*t)
  //          duct = 77.6*(p/t)+373000.*rh*e(t)/(t*t)
  //...................................................................

  if (fieldcalc::bad_hlevel(alevel, blevel)) {
    return false;
  }

  const int fsize = nx * ny;
  const bool inAllDefined = fDefined == miutil::ALL_DEFINED;
  size_t n_undefined = 0;
  MIUTIL_OPENMP_PARALLEL(fsize, for reduction(+:n_undefined))
  for (int i = 0; i < fsize; i++) {
    if (fieldcalc::is_defined(inAllDefined, t[i], h[i], ps[i], undef)) {
      const float p = fieldcalc::p_hlevel(ps[i], alevel, blevel);
      float tk = t[i];
      if (compute % 2 == 0)
        tk *= fieldcalc::pidcp_from_p(p);
      if (compute == 1 or compute == 2) { // T(Kelvin),q -> ducting  or  TH,q -> ducting
        duct[i] = fieldcalc::tk_q_duct(tk, h[i], p);
      } else if (compute == 3 or compute == 4) { // T(Kelvin),RH(%) -> ducting  or TH,RH(%) -> ducting
        duct[i] = fieldcalc::tk_rh_duct(tk, h[i], p, undef, n_undefined);
      }
    } else {
      duct[i] = undef;
      n_undefined += 1;
    }
  }
  fDefined = miutil::checkDefined(n_undefined, fsize);
  return true;
}

bool hlevelpressure(int nx, int ny, const float* ps, float alevel, float blevel, float* p, ValuesDefined& fDefined, float undef)
{
  //  Model levels, eta(hybrid) or norlam_sigma:
  //    eta:   alevel,blevel    p = alevel + blevel * ps[]
  //    sigma: alevel=ptop*(1.-sigma)   blevel=sigma  (p=pt+(ps[]-pt)*sigma)
  //
  //  input:
  //    ps[nx*ny]  : surface pressure (hPa)
  //  output:
  //    p[nx*ny]: pressure (hPa)
  if (fieldcalc::bad_hlevel(alevel, blevel))
    return false;

  const int fsize = nx * ny;

  const bool inAllDefined = fDefined == miutil::ALL_DEFINED;
  size_t n_undefined = 0;
  MIUTIL_OPENMP_PARALLEL(fsize, for reduction(+:n_undefined))
  for (int i = 0; i < fsize; i++) {
    if (fieldcalc::is_defined(inAllDefined, ps[i], undef))
      p[i] = fieldcalc::p_hlevel(ps[i], alevel, blevel);
    else {
      p[i] = undef;
      n_undefined += 1;
    }
  }
  fDefined = miutil::checkDefined(n_undefined, fsize);
  return true;
}

//---------------------------------------------------
// atmospheric model level (ALEVEL) functions
//---------------------------------------------------

bool aleveltemp(int nx, int ny, const float* tinp, const float* p, const std::string& unit, int compute, float* tout, ValuesDefined& fDefined, float undef)
{
  //  Ukjente modell-flater, gitt trykk (p):
  //
  //     compute=1 ; potensiell temp. -> temp. (grader Celsius)
  //     compute=2 ; potensiell temp. -> temp. (grader Kelvin)
  //     compute=3 ; temp. (grader Kelvin)  -> potensiell temp. (Kelvin)
  //     compute=4 ; temp. (grader Kelvin) -> saturated equivalent pot.temp.
  //     compute=5 ; potensiell temp.      -> saturated equivalent pot.temp.
  if (compute <= 0 or compute >= 6)
    return false;

  if (compute < 3) {
    if (unit == "celsius")
      compute = 1;
    else if (unit == "kelvin")
      compute = 2;
  }

  const int fsize = nx * ny;
  const bool inAllDefined = fDefined == miutil::ALL_DEFINED;
  size_t n_undefined = 0;
  MIUTIL_OPENMP_PARALLEL(fsize, for reduction(+:n_undefined))
  for (int i = 0; i < fsize; i++) {
    if (fieldcalc::is_defined(inAllDefined, tinp[i], p[i], undef)) {
      if (compute == 1) { // TH -> T(Celsius)
        tout[i] = tinp[i] * fieldcalc::pidcp_from_p(p[i]) - t0;
      } else if (compute == 2) { // TH -> T(Kelvin)
        tout[i] = tinp[i] * fieldcalc::pidcp_from_p(p[i]);
      } else if (compute == 3) { // T(Kelvin) -> TH
        tout[i] = tinp[i] / fieldcalc::pidcp_from_p(p[i]);
      } else if (compute == 4) { // T(Kelvin) -> THESAT
        tout[i] = fieldcalc::t_thesat(tinp[i], p[i], fieldcalc::pi_from_p(p[i]), undef, n_undefined);
      } else if (compute == 5) { // TH -> THESAT
        tout[i] = fieldcalc::th_thesat(tinp[i], p[i], fieldcalc::pi_from_p(p[i]), undef, n_undefined);
      }
    } else {
      tout[i] = undef;
      n_undefined += 1;
    }
  }
  fDefined = miutil::checkDefined(n_undefined, fsize);
  return true;
}

bool alevelthe(int nx, int ny, const float* t, const float* q, const float* p, int compute, float* the, ValuesDefined& fDefined, float undef)
{
  //  Ukjente modell-flater, gitt trykk (p):
  //
  //    compute=1 ; temp. (Kelvin)  og q     -> THE, ekvivalent pot.temp (Kelvin)
  //    compute=2 ; pot. temp. (Kelvin) og q -> THE, ekvivalent pot.temp (Kelvin)
  //
  //  input:
  //    nx,ny      : field dimensions
  //    t[nx*ny]   : temp. (Kelvin)  or  pot.temp. (Kelvin)
  //    q[nx*ny]   : specific humidity (kg/kg)
  //    p[nx*ny]   : pressure (hPa)
  //  output:
  //    the[nx*ny] : equivalent pot. temp. (Kelvin)

  if (compute != 1 and compute != 2)
    return false;

  const int fsize = nx * ny;

  const bool inAllDefined = fDefined == miutil::ALL_DEFINED;
  size_t n_undefined = 0;
  MIUTIL_OPENMP_PARALLEL(fsize, for reduction(+:n_undefined))
  for (int i = 0; i < fsize; i++) {
    if (fieldcalc::is_defined(inAllDefined, t[i], q[i], p[i], undef)) {
      const float pi = fieldcalc::pi_from_p(p[i]);
      if (compute == 1) // T(Kelvin), q -> THE
        the[i] = (t[i] * cp + q[i] * xlh) / pi;
      else if (compute == 2) // TH, q -> THE
        the[i] = t[i] + q[i] * xlh / pi;
    } else {
      the[i] = undef;
      n_undefined += 1;
    }
  }
  fDefined = miutil::checkDefined(n_undefined, fsize);
  return true;
}

bool alevelhum(int nx, int ny, const float* t, const float* huminp, const float* p,
               const std::string& unit, int compute, float* humout, ValuesDefined& fDefined, float undef)
{
  //  Ukjente modell-flater, gitt trykk (p):
  //
  //     compute=1  : temp. (Kelvin) og spes. fukt. -> rel. fuktighet (%)
  //     compute=2  : pot. temp. og spesifikk fukt. -> rel. fuktighet (%)
  //     compute=3  : temp. (Kelvin) og rel. fukt.  -> spesifikk fukt. (kg/kg)
  //     compute=4  : pot. temp. og  relativ fukt.  -> spesifikk fukt. (kg/kg)
  //     compute=5  : temp. (Kelvin) og spes. fukt. -> duggpunkt, Td (Celsius)
  //     compute=6  : pot. temp. og spesifikk fukt. -> duggpunkt, Td (Celsius)
  //     compute=7  : temp. (Kelvin) og rel. fukt.  -> duggpunkt, Td (Celsius)
  //     compute=8  : pot. temp. og  relativ fukt.  -> duggpunkt, Td (Celsius)
  //     compute=9  : temp. (Kelvin) og spes. fukt. -> duggpunkt, Td (Kelvin)
  //     compute=10 : pot. temp. og spesifikk fukt. -> duggpunkt, Td (Kelvin)
  //     compute=11 : temp. (Kelvin) og rel. fukt.  -> duggpunkt, Td (Kelvin)
  //     compute=12 : pot. temp. og  relativ fukt.  -> duggpunkt, Td (Kelvin)
  //
  //  Note: compute=7,11 are pressure independant (and not needed here)

  if (compute <= 0 or compute >= 13)
    return false;

  if (compute > 8 && unit == "celsius")
    compute -= 4;
  else if (compute > 4 && compute <= 8 && unit == "kelvin")
    compute += 4;

  const int fsize = nx * ny;
  const float tdconv = (compute >= 9) ? t0 : 0;
  const bool inAllDefined = fDefined == miutil::ALL_DEFINED;

  size_t n_undefined = 0;
  MIUTIL_OPENMP_PARALLEL(fsize, for reduction(+:n_undefined))
  for (int i = 0; i < fsize; i++) {
    if (fieldcalc::is_defined(inAllDefined, t[i], huminp[i], undef) and ((compute != 7 and compute != 11) or inAllDefined or p[i] != undef)) {
      if (compute == 1) { // T(Kelvin),q -> RH(%)
        humout[i] = fieldcalc::tk_q_rh(t[i], huminp[i], p[i], undef, n_undefined);
      } else if (compute == 2) { // TH,q -> RH(%)
        const float tk = t[i] * fieldcalc::pidcp_from_p(p[i]);
        humout[i] = fieldcalc::tk_q_rh(tk, huminp[i], p[i], undef, n_undefined);
      } else if (compute == 3) { // T(Kelvin),RH(%) -> q
        humout[i] = fieldcalc::tk_rh_q(t[i], huminp[i], p[i], undef, n_undefined);
      } else if (compute == 4) { // TH,RH(%) -> q
        const float tk = t[i] * fieldcalc::pidcp_from_p(p[i]);
        humout[i] = fieldcalc::tk_rh_q(tk, huminp[i], p[i], undef, n_undefined);
      } else if (compute == 5 or compute == 9) { // T(Kelvin),q -> Td(Celsius)  or  T(Kelvin),q -> Td(Kelvin)
        humout[i] = fieldcalc::tk_q_td(t[i], huminp[i], p[i], tdconv, undef, n_undefined);
      } else if (compute == 6 or compute == 10) { // TH,q -> Td(Celsius)  or  TH,q -> Td(Kelvin)
        const float tk = t[i] * fieldcalc::pidcp_from_p(p[i]);
        humout[i] = fieldcalc::tk_q_td(tk, huminp[i], p[i], tdconv, undef, n_undefined);
      } else if (compute == 7 or compute == 11) { // T(Kelvin),RH(%) -> Td(Celsius)  or  T(Kelvin),RH(%) -> Td(Kelvin)
        humout[i] = fieldcalc::tk_rh_td(t[i], huminp[i], tdconv, undef, n_undefined);
      } else if (compute == 8 or compute == 12) { // TH,RH(%) -> Td(Celsius)  or  TH,RH(%) -> Td(Kelvin)
        const float tk = t[i] * fieldcalc::pidcp_from_p(p[i]);
        humout[i] = fieldcalc::tk_rh_td(tk, huminp[i], tdconv, undef, n_undefined);
      }
    } else {
      humout[i] = undef;
      n_undefined += 1;
    }
  }
  fDefined = miutil::checkDefined(n_undefined, fsize);
  return true;
}

bool alevelducting(int nx, int ny, const float* t, const float* h, const float* p, int compute, float* duct, ValuesDefined& fDefined, float undef)
{
  //  Unknown atmospheric level, with input pressure field:
  //
  //    compute=1 ; temp. (Kelvin)      and q (kg/kg)  -> ducting
  //    compute=2 ; pot. temp. (Kelvin) and q (kg/kg)  -> ducting
  //    compute=3 ; temp. (Kelvin)      and RH(%)      -> ducting
  //    compute=4 ; pot. temp. (Kelvin) and RH(%)      -> ducting
  //
  //  input:
  //    nx,ny      : field dimensions
  //    t[nx*ny]   : temp. (Kelvin)  or  pot.temp. (Kelvin)
  //    h[nx*ny]   : specific humidity (kg/kg) or rel.humidity (%)
  //    p[nx*ny]   : pressure (hPa)
  //  output:
  //    duct[nx*ny]: ducting
  //...................................................................
  //       t in unit Kelvin, p in unit hPa, q in unit kg/kg
  //       duct=77.6*(p/t)+373000.*(q*p)/(eps*t*t)
  //       q*p/eps = rh*qsat*p/eps = rh*(eps*e(t)/p)*p/eps
  //               = rh*e(t) = (e(td)/e(t))*e(t) = e(td)
  //       => duct = 77.6*(p/t)+373000.*e(td)/(t*t)
  //          duct = 77.6*(p/t)+373000.*rh*e(t)/(t*t)
  //...................................................................

  const int fsize = nx * ny;
  const bool inAllDefined = fDefined == miutil::ALL_DEFINED;

  size_t n_undefined = 0;
  MIUTIL_OPENMP_PARALLEL(fsize, for reduction(+:n_undefined))
  for (int i = 0; i < fsize; i++) {
    if (fieldcalc::is_defined(inAllDefined, t[i], h[i], p[i], undef)) {
      float tk = t[i];
      if (compute % 2 == 0)
        tk *= fieldcalc::pidcp_from_p(p[i]);
      if (compute == 1 or compute == 2) { // T(Kelvin),q -> ducting  or  TH,q -> ducting
        duct[i] = fieldcalc::tk_q_duct(tk, h[i], p[i]);
      } else if (compute == 3 or compute == 4) { // T(Kelvin),RH(%) -> ducting  or TH,RH(%) -> ducting
        duct[i] = fieldcalc::tk_rh_duct(tk, h[i], p[i], undef, n_undefined);
      }
    } else {
      duct[i] = undef;
    }
  }
  return true;
}

//---------------------------------------------------
// isentropic level (ILEVEL) function
//---------------------------------------------------

bool ilevelgwind(int nx, int ny, const float* mpot, const float* xmapr, const float* ymapr, const float* fcoriolis, float* ug, float* vg,
                 ValuesDefined& fDefined, float undef)
{
  //  Geostophic wind in isentropic level
  //  (centered differences)
  //
  //  isentropic level, z=Montgomery potential (M=cp*T+g*z)
  //
  //  input:  mpot[nx*ny] : Montgomery potential in isentropic level
  //          xmapr[nx*ny],ymapr[nx*ny] : map ratios (xm/(hx*2),ym/(hy*2))
  //          fcoriolis[nx*ny] : coriolis parameter
  //  output: ug[nx*ny] : x component of geostrophic wind (m/s)
  //          vg[nx*ny] : y component of geostrophic wind (m/s)

  if (nx < 3 || ny < 3)
    return false;
  const int fsize = nx * ny;

  // loop extended, reset bad computations at boundaries later
  const bool inAllDefined = fDefined == miutil::ALL_DEFINED;
  size_t n_undefined = 0;
  MIUTIL_OPENMP_PARALLEL(fsize, for reduction(+:n_undefined))
  for (int i = nx; i < fsize - nx; i++) {
    if (fieldcalc::is_defined(inAllDefined, mpot[i - nx], mpot[i - 1], mpot[i + 1], mpot[i + nx], undef)) {
      ug[i] = -0.5 * ymapr[i] * (mpot[i + nx] - mpot[i - nx]) / fcoriolis[i];
      vg[i] = 0.5 * xmapr[i] * (mpot[i + 1] - mpot[i - 1]) / fcoriolis[i];
    } else {
      ug[i] = undef;
      vg[i] = undef;
      n_undefined += 1;
    }
  }
  fDefined = miutil::checkDefined(n_undefined, fsize);

  // fill in edge values not computed (or badly computed) above
  fillEdges(nx, ny, ug);
  fillEdges(nx, ny, vg);
  return true;
}

//---------------------------------------------------
// ocean depth level (OZLEVEL) functions
//---------------------------------------------------

bool seaSoundSpeed(int nx, int ny, const float* t, const float* s, float z_, int compute, float* soundspeed, ValuesDefined& fDefined, float undef)
{
  //-------------------------------------------------------------
  //       D. Ross: "Revised simplified formulae
  //                 for calculating the speed of
  //                 sound in sea water"
  //       SACLANTCEN Memorandum SM-107 15 Mar 78
  //-------------------------------------------------------------
  //
  //    compute=1 ; temp. (Celsius) and salinity (ppt)
  //    compute=2 ; temp. (Kelvin)  and salinity (ppt)
  //
  //  input:
  //    nx,ny      : field dimensions
  //    t[nx*ny]   : temp. (Kelvin)  or  temp. (Celsius)
  //    s[nx*ny]   : salinity (ppt)
  //    z          : depth in meter (constant)
  //  output:
  //    soundspeed[nx*ny]: sea water sound speed in m/s

  if (compute != 1 and compute != 2)
    return false;

  const int fsize = nx * ny;

  const float tconv = (compute == 1) ? 0 : t0;
  const double Z = fabsf(z_);
  const double Cz = 0.01635 * Z + 0.000000175 * Z * Z;

  const bool inAllDefined = fDefined == miutil::ALL_DEFINED;
  size_t n_undefined = 0;
  MIUTIL_OPENMP_PARALLEL(fsize, for reduction(+:n_undefined))
  for (int i = 0; i < fsize; i++) {
    if (fieldcalc::is_defined(inAllDefined, t[i], s[i], undef)) {
      const float T = t[i] - tconv;
      const float S = s[i];
      const double Ct = 4.565 * T - 0.0517 * T * T + 0.000221 * T * T * T;
      const double Cs = (1.338 - 0.013 * T + 0.0001 * T * T) * (S - 35.0);
      const double speed = 1449.1 + Ct + Cs + Cz;
      soundspeed[i] = float(speed);
    } else {
      soundspeed[i] = undef;
      n_undefined += 1;
    }
  }
  fDefined = miutil::checkDefined(n_undefined, fsize);
  return true;
}

//---------------------------------------------------
// level (vertical coordinate) independant functions
//---------------------------------------------------

bool cvtemp(int nx, int ny, const float* tinp, int compute, float* tout, ValuesDefined& fDefined, float undef)
{
  // compute=1 : Temperature from degrees Kelvin to degrees Celsius
  // compute=2 : Temperature from degrees Celsius to degrees Kelvin
  // compute=3 : Temperature from degrees Kelvin to degrees Celsius
  //             if input seams to be Kelvin
  // compute=4 : Temperature from degrees Celsius to degrees Kelvin
  //             if input seams to be Celsius

  int fsize = nx * ny;
  float tconvert;

  switch (compute) {
  case 1:
    tconvert = -t0;
    break;
  case 2:
    tconvert = +t0;
    break;
  case 3:
    tconvert = -t0;
    break;
  case 4:
    tconvert = +t0;
    break;
  default:
    return false;
  }

  const bool inAllDefined = fDefined == miutil::ALL_DEFINED;

  if (compute == 3 || compute == 4) {
    float tavg = 0.;
    int navg = 0;
    MIUTIL_OPENMP_PARALLEL(fsize, for reduction(+: tavg, navg))
    for (int i = 0; i < fsize; i++) {
      if (fieldcalc::is_defined(inAllDefined, tinp[i], undef)) {
        tavg += tinp[i];
        navg += 1;
      }
    }
    if (navg > 0)
      tavg /= float(navg);

    if ((compute == 3 && tavg < t0 / 2.) || (compute == 4 && tavg > t0 / 2.)) {
      if (&tout[0] != &tinp[0]) {
        MIUTIL_OPENMP_PARALLEL(fsize, for)
        for (int i = 0; i < fsize; i++)
          tout[i] = tinp[i];
      }
      return true;
    }
  }

  size_t n_undefined = 0;
  MIUTIL_OPENMP_PARALLEL(fsize, for reduction(+:n_undefined))
  for (int i = 0; i < fsize; i++) {
    if (fieldcalc::is_defined(inAllDefined, tinp[i], undef))
      tout[i] = tinp[i] + tconvert;
    else {
      tout[i] = undef;
      n_undefined += 1;
    }
  }
  fDefined = miutil::checkDefined(n_undefined, fsize);
  return true;
}

bool abshum(int nx, int ny, const float* t, const float* rhum, float* abshumout, ValuesDefined& fDefined, float undef)
{
  // Equations from Vaisala publication:
  // v = 1 - (T/Tc);
  // ln(Pws/Pc) = (Tc/T) * (C1*v + C2*v^1.5 + C3*v^3 + C4*v^3.5 + C5*v^4 + C6*v^7.5);
  // Pw = Pws * relhum/100;
  // abshum = C * Pw/T;

  // T   = Temperature (Kelvin)
  // Pw  = Vapour pressure (Pa)
  // Pws = Saturation vapor pressure (hPa)

  // Example: T = 20 C, RH = 80%
  // Pw = Pws(20°C) * 80/100 = 18.7 hPa
  // A = 2.16679 * 1870/(273.16 + 20) = 13.82

  // EK Oct. 2012

  // function abshum_out = calc_abshum(in_airtmp,in_relhum)

  // Set constants
  const float C = 2.16679; //(gK/J)
  const float C1 = -7.85951783;
  const float C2 = 1.84408259;
  const float C3 = -11.7866497;
  const float C4 = 22.6807411;
  const float C5 = -15.9618719;
  const float C6 = 1.80122502;
  const float Tc = 647.096; // Critical temperature, Kelvin
  const float Pc = 220640;  // Critical pressure, hPa

  const bool inAllDefined = fDefined == miutil::ALL_DEFINED;
  const int fsize = nx * ny;

  size_t n_undefined = 0;
  MIUTIL_OPENMP_PARALLEL(fsize, for reduction(+:n_undefined))
  for (int i = 0; i < fsize; i++) {
    if (fieldcalc::is_defined(inAllDefined, t[i], rhum[i], undef)) {
      // v = 1 - (T/Tc);
      const float v = 1 - t[i] / Tc, tii = 1 / t[i];

      // ln(Pws/Pc) = (Tc/T) * (C1*v + C2*v^1.5 + C3*v^3 + C4*v^3.5 + C5*v^4 + C6*v^7.5);
      // Solve for Pws:
      // Pws/Pc = exp[(Tc/T) * (C1*v + C2*v^1.5 + C3*v^3 + C4*v^3.5 + C5*v^4 + C6*v^7.5)]
      // Pws = Pc * exp[(Tc/T) * (C1*v + C2*v^1.5 + C3*v^3 + C4*v^3.5 + C5*v^4 + C6*v^7.5)]
      const float v2 = v * v, v3 = v * v2, v4 = v2 * v2, v1_5 = v * sqrt(v), v3_5 = v2 * v1_5, v7_5 = v4 * v3_5;
      float Pws = Pc * exp(Tc * tii * (C1 * v + C2 * v1_5 + C3 * v3 + C4 * v3_5 + C5 * v4 + C6 * v7_5));

      // Pw = Pws * relhum/100;
      float Pw = Pws * rhum[i];

      // abshum = C * Pw/T;
      abshumout[i] = C * Pw * 100 * tii;
    } else {
      abshumout[i] = undef;
      n_undefined += 1;
    }
  }
  fDefined = miutil::checkDefined(n_undefined, fsize);
  return true;
}

bool cvhum(int nx, int ny, const float* t, const float* huminp, const std::string& unit, int compute, float* humout, ValuesDefined& fDefined, float undef)
{
  //     compute=1 : temp. (Kelvin)  og rel. fukt.  -> duggpunkt, Td (Kelvin)
  //     compute=2 : temp. (Kelvin)  og rel. fukt.  -> duggpunkt, Td (Celsius)
  //     compute=3 : temp. (Celsius) og rel. fukt.  -> duggpunkt, Td (Celsius)
  //     compute=4 : temp. og duggpunkt (Kelvin)    -> rel.fuktighet (%)
  //     compute=5 : temp. og duggpunkt (Celsius)   -> rel.fuktighet (%)

  float unit_scale = 100;
  if (compute == 1 && unit == "celsius")
    compute = 2;
  if ((compute == 4 || compute == 5) && unit == "1")
    unit_scale = 1;

  const int fsize = nx * ny;
  const float tconv = (compute == 1 || compute == 2 || compute == 4) ? t0 : 0;
  const float tdconv = (compute == 1) ? t0 : 0;
  const bool inAllDefined = fDefined == miutil::ALL_DEFINED;

  switch (compute) {

  case 1: // T(Kelvin),RH(%)  -> Td(Kelvin)
  case 2: // T(Kelvin),RH(%)  -> Td(Celsius)
  case 3: // T(Celsius),RH(%) -> Td(Celsius)
  {
    size_t n_undefined = 0;
    MIUTIL_OPENMP_PARALLEL(fsize, for reduction(+:n_undefined))
    for (int i = 0; i < fsize; i++) {
      if (fieldcalc::is_defined(inAllDefined, t[i], huminp[i], undef)) {
        const ewt_calculator ewt(t[i] - tconv);
        if (not ewt.defined()) {
          humout[i] = undef;
          n_undefined += 1;
        } else {
          const float et = ewt.value();
          const float rh = fieldcalc::clamp_rh(0.01 * huminp[i]);
          const float etd = rh * et;
          const float tdc = ewt.inverse(etd);
          humout[i] = tdc + tdconv;
        }
      } else {
        humout[i] = undef;
        n_undefined += 1;
      }
    }
    fDefined = miutil::checkDefined(n_undefined, fsize);
    break;
  }

  case 4: // T(Kelvin),Td(Kelvin)   -> RH(%)
  case 5: // T(Celsius),Td(Celsius) -> RH(%)
  {
    size_t n_undefined = 0;
    MIUTIL_OPENMP_PARALLEL(fsize, for reduction(+:n_undefined))
    for (int i = 0; i < fsize; i++) {
      if (fieldcalc::is_defined(inAllDefined, t[i], huminp[i], undef)) {
        const ewt_calculator ewt(t[i] - tconv), ewt2(huminp[i] - tconv);
        if (not(ewt.defined() and ewt2.defined())) {
          humout[i] = undef;
          n_undefined = false;
        } else {
          const float et = ewt.value();
          const float etd = ewt2.value();
          const float rh = etd / et;
          humout[i] = rh * unit_scale;
        }
      } else {
        humout[i] = undef;
        n_undefined += 1;
      }
    }
    fDefined = miutil::checkDefined(n_undefined, fsize);
    break;
  }

  default:
    return false;
  }
  return true;
}

bool vectorabs(int nx, int ny, const float* u, const float* v, float* ff, ValuesDefined& fDefined, float undef)
{
  //  ff = sqrt(u*u+v*v)
  //
  //  input:  u[nx*ny],v[nx*ny] : vector/wind components
  //  output: ff[nx*ny]

  const int fsize = nx * ny;

  const bool inAllDefined = fDefined == miutil::ALL_DEFINED;
  size_t n_undefined = 0;
  MIUTIL_OPENMP_PARALLEL(fsize, for reduction(+:n_undefined))
  for (int i = 0; i < fsize; i++) {
    if (fieldcalc::is_defined(inAllDefined, u[i], v[i], undef))
      ff[i] = miutil::absval(u[i], v[i]);
    else {
      ff[i] = undef;
      n_undefined += 1;
    }
  }
  fDefined = miutil::checkDefined(n_undefined, fsize);
  return true;
}

bool relvort(int nx, int ny, const float* u, const float* v, const float* xmapr, const float* ymapr, float* rvort, ValuesDefined& fDefined, float undef)
{
  //  Relative vorticity from u and v (centered differences)
  //
  //  input:  u[nx*ny],v[nx*ny] : wind components
  //          xmapr[nx*ny],ymapr[nx*ny] : map ratios (xm/(hx*2),ym/(hy*2))
  //  output: rvort[nx*ny]

  if (nx < 3 || ny < 3)
    return false;

  const int fsize = nx * ny;

  // loop extended, reset bad computations at boundaries later
  const bool inAllDefined = fDefined == miutil::ALL_DEFINED;
  size_t n_undefined = 0;
  MIUTIL_OPENMP_PARALLEL(fsize-2*nx, for reduction(+:n_undefined))
  for (int i = nx; i < fsize - nx; i++) {
    if (fieldcalc::is_defined(inAllDefined, v[i - 1], v[i + 1], u[i - nx], u[i + nx], undef))
      rvort[i] = 0.5 * xmapr[i] * (v[i + 1] - v[i - 1]) - 0.5 * ymapr[i] * (u[i + nx] - u[i - nx]);
    else {
      rvort[i] = undef;
      n_undefined += 1;
    }
  }
  fDefined = miutil::checkDefined(n_undefined, fsize - 2 * nx);

  // fill in edge values not computed (or badly computed) above
  fillEdges(nx, ny, rvort);
  return true;
}

bool absvort(int nx, int ny, const float* u, const float* v, const float* xmapr, const float* ymapr, const float* fcoriolis, float* avort,
             ValuesDefined& fDefined, float undef)
{
  //  Absolute vorticity from u and v (centered differences)
  //
  //  input:  u[nx*ny],v[nx*ny] : wind components
  //          xmapr[nx*ny],ymapr[nx*ny] : map ratios (xm/(hx*2),ym/(hy*2))
  //          fcoriolis[nx*ny] : coriolis parameter
  //  output: avort[nx*ny]

  if (nx < 3 || ny < 3)
    return false;

  const int fsize = nx * ny;

  // loop extended, reset bad computations at boundaries later
  const bool inAllDefined = fDefined == miutil::ALL_DEFINED;
  size_t n_undefined = 0;
  MIUTIL_OPENMP_PARALLEL(fsize-2*nx, for reduction(+:n_undefined))
  for (int i = nx; i < fsize - nx; i++) {
    if (fieldcalc::is_defined(inAllDefined, v[i - 1], v[i + 1], u[i - nx], u[i + nx], undef))
      avort[i] = 0.5 * xmapr[i] * (v[i + 1] - v[i - 1]) - 0.5 * ymapr[i] * (u[i + nx] - u[i - nx]) + fcoriolis[i];
    else {
      avort[i] = undef;
      n_undefined += 1;
    }
  }
  fDefined = miutil::checkDefined(n_undefined, fsize - 2 * nx);

  // fill in edge values not computed (or badly computed) above
  fillEdges(nx, ny, avort);

  return true;
}

bool divergence(int nx, int ny, const float* u, const float* v, const float* xmapr, const float* ymapr, float* diverg, ValuesDefined& fDefined, float undef)
{
  //  Divergence from u and v (centered differences)
  //
  //  input:  u[nx*ny],v[nx*ny] : wind components
  //          xmapr[nx*ny],ymapr[nx*ny] : map ratios (xm/(hx*2),ym/(hy*2))
  //  output: diverg[nx*ny]

  if (nx < 3 || ny < 3)
    return false;

  const int fsize = nx * ny;
  // loop extended, reset bad computations at boundaries later
  const bool inAllDefined = fDefined == miutil::ALL_DEFINED;
  size_t n_undefined = 0;
  MIUTIL_OPENMP_PARALLEL(fsize-2*nx, for reduction(+:n_undefined))
  for (int i = nx; i < fsize - nx; i++) {
    if (fieldcalc::is_defined(inAllDefined, v[i - 1], v[i + 1], u[i - nx], u[i + nx], undef))
      diverg[i] = 0.5 * xmapr[i] * (u[i + 1] - u[i - 1]) + 0.5 * ymapr[i] * (v[i + nx] - v[i - nx]);
    else {
      diverg[i] = undef;
      n_undefined += 1;
    }
  }
  fDefined = miutil::checkDefined(n_undefined, fsize - 2 * nx);

  // fill in edge values not computed (or badly computed) above
  fillEdges(nx, ny, diverg);

  return true;
}

bool advection(int nx, int ny, const float* f, const float* u, const float* v, const float* xmapr, const float* ymapr, float hours, float* advec,
               ValuesDefined& fDefined, float undef)
{
  // compute advection of a scalar field f (e.g. temperature):
  //     adection = -(u*dt/dx+v*dt/dy)*hours
  //
  //     dt/dx and dt/dy as centered differences
  //
  // input: f[nx*ny]: scalar field (e.g. temperature)
  //        u[nx*ny]: wind component in x-direction (m/s)
  //        v[nx*ny]: wind component in y-direction (m/s)
  //        xmapr[nx*ny],ymapr[nx*ny] : map ratios (xm/(hx*2),ym/(hy*2))
  //        hours: scaling to get a chosen time unit (delta_f/timeunit)
  //               e.g. hours=1./3600. => per second. ('SI')
  //                    hours=1.       => per hour
  //                    hours=24.      => per 24 hours
  // output: advec[nx*ny]

  if (nx < 3 || ny < 3)
    return false;

  const float scale = -3600. * hours;
  const int fsize = nx * ny;

  // loop extended, reset bad computations at boundaries later
  const bool inAllDefined = fDefined == miutil::ALL_DEFINED;
  size_t n_undefined = 0;
  MIUTIL_OPENMP_PARALLEL(fsize-2*nx, for reduction(+:n_undefined))
  for (int i = nx; i < fsize - nx; i++) {
    if (fieldcalc::is_defined(inAllDefined, u[i], v[i], f[i - nx], f[i - 1], f[i + 1], f[i + nx], undef))
      advec[i] = (u[i] * 0.5 * xmapr[i] * (f[i + 1] - f[i - 1]) + v[i] * 0.5 * ymapr[i] * (f[i + nx] - f[i - nx])) * scale;
    else {
      advec[i] = undef;
      n_undefined += 1;
    }
  }
  fDefined = miutil::checkDefined(n_undefined, fsize - 2 * nx);

  // fill in edge values not computed (or badly computed) above
  fillEdges(nx, ny, advec);
  return true;
}

bool gradient(int nx, int ny, const float* field, const float* xmapr, const float* ymapr, int compute, float* fgrad, ValuesDefined& fDefined, float undef)
{
  //  Gradienter beregnet med sentrerte differanser, kartfaktor
  //  og gridavstand i meter.
  //
  //  compute=1: df/dx
  //  compute=2: df/dy
  //  compute=3: abs(del(f))= sqrt((df/dx)**2 + (df/dy)**2)
  //  compute=4: del2(f)= del(del(f))
  //
  //  xmapr= xm/(hx*2)
  //  ymapr= ym/(hy*2)
  //
  //  NB! for compute=4 (del2(f)) behandles kartfaktorene ikke korrekt
  //      siden den deriverte i x- og y-retning burde inngaa,
  //      komplisert (mye regning) og antas aa ha liten effekt !!!

  const int fsize = nx * ny;

  if (nx < 3 || ny < 3)
    return false;

  const bool inAllDefined = fDefined == miutil::ALL_DEFINED;
  size_t n_undefined = 0;

  switch (compute) {

  case 1: // df/dx (loop extended, reset/complete later
    MIUTIL_OPENMP_PARALLEL(fsize-2, for reduction(+:n_undefined))
    for (int i = 1; i < fsize - 1; i++) {
      if (fieldcalc::is_defined(inAllDefined, field[i - 1], field[i + 1], undef))
        fgrad[i] = 0.5 * xmapr[i] * (field[i + 1] - field[i - 1]);
      else {
        fgrad[i] = undef;
        n_undefined += 1;
      }
    }
    break;

  case 2: // df/dy (loop extended, reset/complete later
    MIUTIL_OPENMP_PARALLEL(fsize-2*nx, for reduction(+:n_undefined))
    for (int i = nx; i < fsize - nx; i++) {
      if (fieldcalc::is_defined(inAllDefined, field[i - nx], field[i + nx], undef))
        fgrad[i] = 0.5 * ymapr[i] * (field[i + nx] - field[i - nx]);
      else {
        fgrad[i] = undef;
        n_undefined += 1;
      }
    }
    break;

  case 3: // abs(del(f))= sqrt((df/dx)**2 + (df/dy)**2)
    MIUTIL_OPENMP_PARALLEL(fsize-2*nx, for reduction(+:n_undefined))
    for (int i = nx; i < fsize - nx; i++) {
      if (fieldcalc::is_defined(inAllDefined, field[i - nx], field[i - 1], field[i + 1], field[i + nx], undef)) {
        const float dfdx = 0.5 * xmapr[i] * (field[i + 1] - field[i - 1]);
        const float dfdy = 0.5 * ymapr[i] * (field[i + nx] - field[i - nx]);
        fgrad[i] = miutil::absval(dfdx, dfdy);
      } else {
        fgrad[i] = undef;
        n_undefined += 1;
      }
    }
    break;

  case 4: // delsquare(f)= del(del(f))
    MIUTIL_OPENMP_PARALLEL(fsize-2*nx, for reduction(+:n_undefined))
    for (int i = nx; i < fsize - nx; i++) {
      if (fieldcalc::is_defined(inAllDefined, field[i - nx], field[i - 1], field[i], field[i + 1], field[i + nx], undef)) {
        const float d2fdx = field[i - 1] - 2.0 * field[i] + field[i + 1];
        const float d2fdy = field[i - nx] - 2.0 * field[i] + field[i + nx];
        fgrad[i] = 4.0 * (0.25 * xmapr[i] * xmapr[i] * d2fdx + 0.25 * ymapr[i] * ymapr[i] * d2fdy);
      } else {
        fgrad[i] = undef;
        n_undefined += 1;
      }
    }
    break;

  default:
    return false;
  }

  fDefined = miutil::checkDefined(n_undefined, fsize - 2 * nx);

  // fill in edge values not computed (or badly computed) above
  fillEdges(nx, ny, fgrad);

  return true;
}

bool shapiro2_filter(int nx, int ny, float* field, float* fsmooth, ValuesDefined& fDefined, float undef)
{
  //  2nd order Shapiro filter
  //
  //  G.J.Haltiner, Numerical Weather Prediction,
  //                   Objective Analysis,
  //                       Smoothing and filtering
  //
  //  input:   nx,ny          - field dimensions
  //           field[nx*ny]   - the field
  //  output:  fsmooth[nx*ny] - the smoothed field
  //
  //  Note: field and fsmooth may be the same array!

  const int fsize = nx * ny;
  float s;

  if (nx < 3 || ny < 3)
    return false;

  float* f1 = 0;
  float* f2 = new float[fsize];

  if (field == fsmooth) {
    f1 = field;
  } else {
    for (int i = 0; i < fsize; i++)
      fsmooth[i] = field[i];
    f1 = fsmooth;
  }

  const bool inAllDefined = fDefined == miutil::ALL_DEFINED;
  if (inAllDefined) {

    s = 0.25;

    for (int n = 0; n < 2; n++) {

      for (int i = 1; i < fsize - 1; i++)
        f2[i] = f1[i] + s * (f1[i - 1] + f1[i + 1] - 2. * f1[i]);

      for (int j = 0, i1 = 0, i2 = nx - 1; j < ny; j++, i1 += nx, i2 += nx) {
        f2[i1] = f1[i1];
        f2[i2] = f1[i2];
      }

      for (int i = nx; i < fsize - nx; i++)
        f1[i] = f2[i] + s * (f2[i - nx] + f2[i + nx] - 2. * f2[i]);

      for (int i1 = 0, i2 = fsize - nx; i1 < nx; i1++, i2++) {
        f1[i1] = f2[i1];
        f1[i2] = f2[i2];
      }

      s = -0.25;
    }

  } else {

    // s is uninitialised at this point!
    s = 0.25; // moved it here...

    float* s1 = new float[fsize];
    float* s2 = new float[fsize];

    for (int i = 1; i < fsize - 1; i++)
      s1[i] = fieldcalc::is_defined(inAllDefined, f1[i - 1], f1[i], f1[i + 1], undef) ? s : 0;

    for (int i = nx; i < fsize - nx; i++)
      s2[i] = fieldcalc::is_defined(inAllDefined, f1[i - nx], f1[i], f1[i + nx], undef) ? s : 0;

    // s = 0.25; // ... from here

    for (int n = 0; n < 2; n++) {

      for (int i = 1; i < fsize - 1; i++)
        f2[i] = f1[i] + s1[i] * (f1[i - 1] + f1[i + 1] - 2 * f1[i]);

      for (int j = 0, i1 = 0, i2 = nx - 1; j < ny; j++, i1 += nx, i2 += nx) {
        f2[i1] = f1[i1];
        f2[i2] = f1[i2];
      }

      for (int i = nx; i < fsize - nx; i++)
        f1[i] = f2[i] + s2[i] * (f2[i - nx] + f2[i + nx] - 2 * f2[i]);

      for (int i1 = 0, i2 = fsize - nx; i1 < nx; i1++, i2++) {
        f1[i1] = f2[i1];
        f1[i2] = f2[i2];
      }

      s = -0.25;
    }

    delete[] s1;
    delete[] s2;
  }

  delete[] f2;

  fDefined = miutil::ALL_DEFINED;

  return true;
}

bool windCooling(int nx, int ny, const float* t, const float* u, const float* v, int compute, float* dtcool, ValuesDefined& fDefined, float undef)
{
  // Compute wind cooling (the difference, not the sensed temperature)
  //
  //  compute=1 : temp. (Kelvin)
  //  compute=2 : temp. (Celsius)
  //
  //  input:   nx,ny         - field dimensions
  //           t[nx*ny]      - temperature (Kelvin or Celsius)
  //           u[nx*ny]      - wind component x-direction (m/s)
  //           v[nx*ny]      - wind component y-direction (m/s)
  //  output:  dtcool[nx*ny] - the cooling (temp. difference)

  int fsize = nx * ny;
  float tconv;

  tconv = 0.;
  if (compute == 1)
    tconv = t0;

  switch (compute) {

  case 1: // T(Kelvin)  -> dTcooling
  case 2: // T(Celsius) -> dTcooling
  {
    const bool inAllDefined = fDefined == miutil::ALL_DEFINED;
    size_t n_undefined = 0;
    MIUTIL_OPENMP_PARALLEL(fsize, for reduction(+:n_undefined))
    for (int i = 0; i < fsize; i++) {
      if (fieldcalc::is_defined(inAllDefined, t[i], u[i], v[i], undef)) {
        const float tc = t[i] - tconv;
        const float ff = miutil::absval(u[i], v[i]) * 3.6; // m/s -> km/h
        const float ffpow = powf(ff, 0.16);
        dtcool[i] = 13.12 + 0.6215 * tc - 11.37 * ffpow + 0.3965 * tc * ffpow;
        if (dtcool[i] > 0.)
          dtcool[i] = 0.;
      } else {
        dtcool[i] = undef;
        n_undefined += 1;
      }
    }
    break;
  }

  default:
    return false;
  }
  return true;
}

bool underCooledRain(int nx, int ny, const float* precip, const float* snow, const float* tk, float precipMin, float snowRateMax,
                     float tcMax, float* undercooled, ValuesDefined& fDefined, float undef)
{
  // compute possibility/danger for undercooled rain.....
  //
  // input:  precip[]  : total precip (rain+snow) in unit mm
  //         snow[]    : snow in unit mm (water)
  //         tk[]      : temperature im unit Kelvin (t0m)
  //         undef     : unde|fined value for all fields
  //         maxSnowRate: max rate of snow in total precip
  //         precLimit : min precipitation in unit mm
  //         tempLimit : max temperature in unit Celsius
  // output: undercooled[] : 0=no 1=yes

  const int fsize = nx * ny;
  const float tkMax = tcMax + t0;

  const bool inAllDefined = fDefined == miutil::ALL_DEFINED;
  size_t n_undefined = 0;
  MIUTIL_OPENMP_PARALLEL(fsize, for reduction(+:n_undefined))
  for (int i = 0; i < fsize; i++) {
    if (fieldcalc::is_defined(inAllDefined, precip[i], snow[i], tk[i], undef)) {
      if (precip[i] >= precipMin && tk[i] <= tkMax && snow[i] <= precip[i] * snowRateMax)
        undercooled[i] = 1.;
      else
        undercooled[i] = 0.;
    } else {
      undercooled[i] = undef;
      n_undefined += 1;
    }
  }
  fDefined = miutil::checkDefined(n_undefined, fsize);
  return true;
}

bool thermalFrontParameter(int nx, int ny, const float* tx, const float* xmapr, const float* ymapr, float* tfp, ValuesDefined& fDefined, float undef)
{
  // Thermal Front Parameter
  // TFP = -del(abs(del(T)) . del(t)/abs(del(T))
  //
  // input:   tx[nx*ny] : temperature, pot.temp. or some other temperature
  // output: tfp[nx*ny]
  //
  //  xmapr= xm/(hx*2)
  //  ymapr= ym/(hy*2)

  const int fsize = nx * ny;

  std::unique_ptr<float[]> absdelt(new float[fsize]);

  if (!gradient(nx, ny, tx, xmapr, ymapr, 3, absdelt.get(), fDefined, undef)) {
    return false;
  }

  // loop extended, reset bad computations at boundaries later
  const bool inAllDefined = fDefined == miutil::ALL_DEFINED;
  size_t n_undefined = 0;
  MIUTIL_OPENMP_PARALLEL(fsize-2*nx, for reduction(+:n_undefined))
  for (int i = nx; i < fsize - nx; i++) {
    if (fieldcalc::is_defined(inAllDefined, tx[i - nx], tx[i - 1], tx[i + 1], tx[i + nx], absdelt[i - nx], absdelt[i - 1], absdelt[i], absdelt[i + 1],
                                 absdelt[i + nx], undef) &&
        absdelt[i] != 0) {
      const float dabsdeltdx = 0.5 * xmapr[i] * (absdelt[i + 1] - absdelt[i - 1]);
      const float dabsdeltdy = 0.5 * ymapr[i] * (absdelt[i + nx] - absdelt[i - nx]);
      const float dtdxa = 0.5 * xmapr[i] * (tx[i + 1] - tx[i - 1]) / absdelt[i];
      const float dtdya = 0.5 * ymapr[i] * (tx[i + nx] - tx[i - nx]) / absdelt[i];
      tfp[i] = -(dabsdeltdx * dtdxa + dabsdeltdy * dtdya);
    } else {
      tfp[i] = undef;
      n_undefined += 1;
    }
  }
  fDefined = miutil::checkDefined(n_undefined, fsize - 2 * nx);

  // fill in edge values not computed (or badly computed) above
  fillEdges(nx, ny, tfp);

  return true;
}

bool pressure2FlightLevel(int nx, int ny, const float* pressure, float* flightlevel, ValuesDefined& fDefined, float undef)
{
  // compute standard FlightLevel from pressure,
  // using the same table as used when mapping standard pressure levels
  // to standard FlightLevels
  // (not very correct, but consistent...)

  // input:  pressure[] in unit hPa
  // output: flightlevel[] in unit 100 feet

  // tables in diMetConstants.h

  const int nTab = nLevelTable - 1;
  const int fsize = nx * ny;

  const bool inAllDefined = fDefined == miutil::ALL_DEFINED;
  size_t n_undefined = 0;
  MIUTIL_OPENMP_PARALLEL(fsize, for reduction(+:n_undefined))
  for (int i = 0; i < fsize; i++) {
    if (fieldcalc::is_defined(inAllDefined, pressure[i], undef)) {
      // TODO use binary search, set k immediately for p outside pLevelTable
      float p = pressure[i];
      if (p > pLevelTable[0])
        p = pLevelTable[0];
      if (p < pLevelTable[nTab])
        p = pLevelTable[nTab];
      int k = 1;
      while (k < nTab && pLevelTable[k] > p)
        k++;
      const float ratio = (p - pLevelTable[k - 1]) / (pLevelTable[k] - pLevelTable[k - 1]);
      flightlevel[i] = fLevelTable[k - 1] + (fLevelTable[k] - fLevelTable[k - 1]) * ratio;
    } else {
      flightlevel[i] = undef;
      n_undefined += 1;
    }
  }
  fDefined = miutil::checkDefined(n_undefined, fsize);
  return true;
}

bool momentumXcoordinate(int nx, int ny, const float* v, const float* xmapr, const float* fcoriolis, float fcoriolisMin, float* mxy, ValuesDefined& fDefined,
                         float undef)
{
  // Momentum coordinates, X component
  // Compute m(x,y) from v (centered differences)
  //
  // fcpriolisMin : minimum value of coriolis parameter used in computation
  //
  // m(x,y) = x + v*xm/(hx*fc) ..... xmapr = xm/(2*hx) !!!

  const int fsize = nx * ny;

  if (nx < 3 || ny < 3)
    return false;

  const float fcormin = fabsf(fcoriolisMin), fcormax = -fcormin;

  const bool inAllDefined = fDefined == miutil::ALL_DEFINED;
  size_t n_undefined = 0;
  MIUTIL_OPENMP_PARALLEL(fsize, for reduction(+:n_undefined))
  for (int i = 0; i < fsize; i++) {
    if (fieldcalc::is_defined(inAllDefined, v[i], undef)) {
      float fcor = fcoriolis[i];
      if (fcor >= 0. && fcor < fcormin)
        fcor = fcormin;
      else if (fcor <= 0. && fcor > fcormax)
        fcor = fcormax;
      mxy[i] = float(i % nx) + v[i] * xmapr[i] / fcor;
    } else {
      mxy[i] = undef;
      n_undefined += 1;
    }
  }
  fDefined = miutil::checkDefined(n_undefined, fsize);
  return true;
}

bool momentumYcoordinate(int nx, int ny, const float* u, const float* ymapr, const float* fcoriolis, float fcoriolisMin, float* nxy, ValuesDefined& fDefined,
                         float undef)
{
  // Momentum coordinates, Y component
  // Compute n(x,y) from u (centered differences)
  //
  // fcpriolisMin : minimum value of coriolis parameter used in computation
  //
  // n(x,y) = y - u*ym/(hy*fc) ..... ymapr = ym/(2*hy) !!!

  if (nx < 3 || ny < 3)
    return false;

  const int fsize = nx * ny;
  const float fcormin = fabsf(fcoriolisMin), fcormax = -fcormin;

  const bool inAllDefined = fDefined == miutil::ALL_DEFINED;
  size_t n_undefined = 0;
  MIUTIL_OPENMP_PARALLEL(fsize, for reduction(+:n_undefined))
  for (int i = 0; i < fsize; i++) {
    if (fieldcalc::is_defined(inAllDefined, u[i], undef)) {
      float fcor = fcoriolis[i];
      if (fcor >= 0. && fcor < fcormin)
        fcor = fcormin;
      else if (fcor <= 0. && fcor > fcormax)
        fcor = fcormax;
      nxy[i] = float(i / nx) - u[i] * ymapr[i] / fcor;
    } else {
      nxy[i] = undef;
      n_undefined += 1;
    }
  }
  fDefined = miutil::checkDefined(n_undefined, fsize);
  return true;
}

bool jacobian(int nx, int ny, const float* field1, const float* field2, const float* xmapr, const float* ymapr, float* fjacobian, ValuesDefined& fDefined,
              float undef)
{
  //  Beregner den jacobiske av f1 og f2:
  //       Jacobian = df1/dx * df2/dy - df1/dy * df2/dx
  //
  //   xmapr = xm/(hx*2)
  //   ymapr = ym/(hy*2)

  if (nx < 3 || ny < 3)
    return false;

  const int fsize = nx * ny;

  // loop extended, reset bad computations at boundaries later
  const bool inAllDefined = fDefined == miutil::ALL_DEFINED;
  size_t n_undefined = 0;
  MIUTIL_OPENMP_PARALLEL(fsize-2*nx, for reduction(+:n_undefined))
  for (int i = nx; i < fsize - nx; i++) {
    if (fieldcalc::is_defined(inAllDefined, field1[i - nx], field1[i - 1], field1[i + 1], field1[i + nx], field2[i - nx], field2[i - 1], field2[i + 1],
                                 field2[i + nx], undef)) {
      const float df1dx = 0.5 * xmapr[i] * (field1[i + 1] - field1[i - 1]);
      const float df1dy = 0.5 * ymapr[i] * (field1[i + nx] - field1[i - nx]);
      const float df2dx = 0.5 * xmapr[i] * (field2[i + 1] - field2[i - 1]);
      const float df2dy = 0.5 * ymapr[i] * (field2[i + nx] - field2[i - nx]);
      fjacobian[i] = df1dx * df2dy - df1dy * df2dx;
    } else {
      fjacobian[i] = undef;
      n_undefined += 1;
    }
  }
  fDefined = miutil::checkDefined(n_undefined, fsize - 2 * nx);

  // fill in edge values not computed (or badly computed) above
  fillEdges(nx, ny, fjacobian);
  return true;
}

bool values2classes(int nx, int ny, const float* fvalue, float* fclass, const std::vector<float>& values, ValuesDefined& fDefined, float undef)
{
  // From value ranges to classes
  // values not covered for will be set to undefined
  //
  // Note: The classes will always be 0,1,2,3,...
  //       The values are minimum values except the last which is a maximim.
  //       If there are 10 values they can become classes 0-8 !
  //       Ranges are [min,max).
  //
  // Note: The argument 'allDefined' may be changed!

  const int fsize = nx * ny;

  if (values.size() < 2)
    return false;

  const int nvalues = values.size() - 2;
  const float fmin = values[0];
  const float fmax = values[nvalues + 1];

  const bool inAllDefined = fDefined == miutil::ALL_DEFINED;
  size_t n_undefined = 0;
  MIUTIL_OPENMP_PARALLEL(fsize, for reduction(+:n_undefined))
  for (int i = 0; i < fsize; i++) {
    if (fieldcalc::is_defined(inAllDefined, fvalue[i], undef) && fvalue[i] >= fmin && fvalue[i] < fmax) {
      int j = 1;
      while (j < nvalues && values[j] < fvalue[i])
        j++;
      fclass[i] = float(j - 1);
    } else {
      fclass[i] = undef;
      n_undefined += 1;
    }
  }
  fDefined = miutil::checkDefined(n_undefined, fsize);
  return true;
}

void minvalueFields(int nx, int ny, const float* field1, const float* field2, float* fres, ValuesDefined& fDefined, float undef)
{
  // use lambda becaus std::min is inline
  binaryFunctionFieldField([](float a, float b) { return std::min(a, b); }, nx, ny, field1, field2, fres, fDefined, undef);
}

void minvalueFieldConst(int nx, int ny, const float* field1, const float value, float* fres, ValuesDefined& fDefined, float undef)
{
  if (value == undef) {
    fillUndef(nx, ny, fres, fDefined, undef);
  } else {
    unaryFunctionField([value](float a) { return std::min(a, value); }, nx, ny, field1, fres, fDefined, undef);
  }
}

void maxvalueFields(int nx, int ny, const float* field1, const float* field2, float* fres, ValuesDefined& fDefined, float undef)
{
  // use lambda becaus std::max is inline
  binaryFunctionFieldField([](float a, float b) { return std::max(a, b); }, nx, ny, field1, field2, fres, fDefined, undef);
}

void maxvalueFieldConst(int nx, int ny, const float* field1, const float value, float* fres, ValuesDefined& fDefined, float undef)
{
  if (value == undef) {
    fillUndef(nx, ny, fres, fDefined, undef);
  } else {
    unaryFunctionField([value](float a) { return std::max(a, value); }, nx, ny, field1, fres, fDefined, undef);
  }
}

void absvalueField(int nx, int ny, const float* field, float* fres, ValuesDefined& fDefined, float undef)
{
  unaryFunctionField(std::abs, nx, ny, field, fres, fDefined, undef);
}

void log10Field(int nx, int ny, const float* field, float* fres, ValuesDefined& fDefined, float undef)
{
  unaryFunctionField(std::log10, nx, ny, field, fres, fDefined, undef);
}

void pow10Field(int nx, int ny, const float* field, float* fres, ValuesDefined& fDefined, float undef)
{
  unaryFunctionField(miutil::pow10, nx, ny, field, fres, fDefined, undef);
}

void logField(int nx, int ny, const float* field, float* fres, ValuesDefined& fDefined, float undef)
{
  unaryFunctionField(std::log, nx, ny, field, fres, fDefined, undef);
}

void expField(int nx, int ny, const float* field, float* fres, ValuesDefined& fDefined, float undef)
{
  unaryFunctionField(std::exp, nx, ny, field, fres, fDefined, undef);
}

void powerField(int nx, int ny, const float* field, float value, float* fres, ValuesDefined& fDefined, float undef)
{
  if (value == undef) {
    fillUndef(nx, ny, fres, fDefined, undef);
  } else {
    unaryFunctionField([value](float a) { return std::pow(a, value); }, nx, ny, field, fres, fDefined, undef);
  }
}

void replaceUndefined(int nx, int ny, const float* field, float value, float* fres, ValuesDefined& fDefined, float undef)
{
  if (value == undef || fDefined == miutil::ALL_DEFINED) {
    // replace all undefined with undef => no change
    copy_field(fres, field, nx*ny);
    return;
  }

  if (fDefined == miutil::NONE_DEFINED) {
    // none defined => all replaced
    std::fill(fres, fres + nx*ny, value);
  } else {
    // SOME_DEFINED
    const size_t fsize = nx * ny;
    MIUTIL_OPENMP_PARALLEL(fsize, for)
    for (size_t i = 0; i < fsize; i++)
      fres[i] = (field[i] == undef) ? value : field[i];
  }
  // ALL_DEFINED as value != undef
  fDefined = miutil::ALL_DEFINED;
}

void replaceDefined(int nx, int ny, const float* field, float value, float* fres, ValuesDefined& fDefined, float undef)
{
  if (value == undef || fDefined == miutil::NONE_DEFINED) {
    // replace all defined with undef, or none defined => will bo NONE_DEFINED afterwards
    std::fill(fres, fres + nx*ny, undef);
    fDefined = miutil::NONE_DEFINED;
    return;
  }

  if (fDefined == miutil::ALL_DEFINED) {
    // all defined => all replaced
    std::fill(fres, fres + nx*ny, value);
  } else {
    // SOME_DEFINED
    const size_t fsize = nx * ny;
    MIUTIL_OPENMP_PARALLEL(fsize, for)
    for (size_t i = 0; i < fsize; i++)
      fres[i] = (field[i] != undef) ? value : field[i];
  }
  // ALL_DEFINED as value != undef
  fDefined = miutil::ALL_DEFINED;
}

//! field1 <+-*/> field2
bool fieldOPERfield(int compute, int nx, int ny, const float* field1, const float* field2, float* fres, ValuesDefined& fDefined, float undef)
{
  switch (compute) {
  case 1: //  compute=1 : field1 + field2
    return binaryFunctionFieldField(std::plus<float>(), nx, ny, field1, field2, fres, fDefined, undef);
  case 2: // field1 - field2
    return binaryFunctionFieldField(std::minus<float>(), nx, ny, field1, field2, fres, fDefined, undef);
  case 3: // field1 * field2
    return binaryFunctionFieldField(std::multiplies<float>(), nx, ny, field1, field2, fres, fDefined, undef);
  case 4: // field1 / field2
    return binaryFunctionFieldFieldUndef(divideUndef, nx, ny, field1, field2, fres, fDefined, undef);
  default:
    return false;
  }
}

bool fieldOPERconstant(int compute, int nx, int ny, const float* field, float value, float* fres, ValuesDefined& fDefined, float undef)
{
  if ((value == undef) || (compute == 4 && value == 0)) // field / constant
    return fillUndef(nx, ny, fres, fDefined, undef);

  // field <+-*/> constant
  switch (compute) {
  case 1: // field + constant
    return unaryFunctionField([value](float f) { return f + value; }, nx, ny, field, fres, fDefined, undef);
  case 2: // field - constant
    return unaryFunctionField([value](float f) { return f - value; }, nx, ny, field, fres, fDefined, undef);
  case 3: // field * constant
    return unaryFunctionField([value](float f) { return f * value; }, nx, ny, field, fres, fDefined, undef);
  case 4: // field / constant
    return unaryFunctionField([value](float f) { return f / value; }, nx, ny, field, fres, fDefined, undef);
  default:
    return false;
  }
}

bool constantOPERfield(int compute, int nx, int ny, float value, const float* field, float* fres, ValuesDefined& fDefined, float undef)
{
  // constant <+-*/> field

  if (value == undef)
    return fillUndef(nx, ny, fres, fDefined, undef);

  // field <+-*/> constant
  switch (compute) {
  case 1: // constant + field
    return unaryFunctionField([value](float f) { return value + f; }, nx, ny, field, fres, fDefined, undef);
  case 2: // constant - field
    return unaryFunctionField([value](float f) { return value - f; }, nx, ny, field, fres, fDefined, undef);
  case 3: // constant * field
    return unaryFunctionField([value](float f) { return value * f; }, nx, ny, field, fres, fDefined, undef);
  case 4: { // constant / field
    auto div = [value](float f, float undef, size_t& n_undefined) { return divideUndef(value, f, undef, n_undefined); };
    return unaryFunctionFieldUndef(div, nx, ny, field, fres, fDefined, undef);
  }
  default:
    return false;
  }
}

bool sumFields(int nx, int ny, const std::vector<float*>& fields, float* fres, ValuesDefined& fDefined, float undef)
{
  // field + field + field +

  const int fsize = nx * ny;
  const size_t nFields = fields.size();
  const bool inAllDefined = fDefined == miutil::ALL_DEFINED;
  size_t n_undefined = 0;
  MIUTIL_OPENMP_PARALLEL(fsize, for reduction(+:n_undefined))
  for (int i = 0; i < fsize; i++) {
    fres[i] = 0;
    for (size_t j = 0; j < nFields; j++) {
      if (fieldcalc::is_defined(inAllDefined, fields[j][i], undef)) {
        fres[i] += fields[j][i];
      } else {
        fres[i] = undef;
        n_undefined += 1;
        break;
      }
    }
  }
  fDefined = miutil::checkDefined(n_undefined, fsize);
  return true;
}

bool meanValue(int nx, int ny, const std::vector<float*>& fields, const std::vector<ValuesDefined>& fDefinedIn, float* fres, ValuesDefined& fDefinedOut,
               float undef)
{
  // returns mean value in each grid point
  // (field + field + field +) / nFields

  const int fsize = nx * ny;
  const size_t nFields = fields.size();
  size_t n_undefined = 0;
  MIUTIL_OPENMP_PARALLEL(fsize, for reduction(+:n_undefined))
  for (int i = 0; i < fsize; i++) {
    fres[i] = 0;
    int nfields_defined = 0;
    for (size_t j = 0; j < nFields; j++) {
      if (fieldcalc::is_defined(fDefinedIn[j] == miutil::ALL_DEFINED, fields[j][i], undef)) {
        nfields_defined++;
        fres[i] += fields[j][i];
      }
    }
    if (nfields_defined > 0) {
      fres[i] /= nfields_defined;
    } else {
      fres[i] = undef;
      n_undefined += 1;
    }
  }
  fDefinedOut = miutil::checkDefined(n_undefined, fsize);
  return true;
}

bool stddevValue(int nx, int ny, const std::vector<float*>& fields, const std::vector<ValuesDefined>& fDefinedIn, float* fres, ValuesDefined& fDefinedOut,
                 float undef)
{
  // returns mean value in each grid point
  // (field + field + field +) / nFields

  const int fsize = nx * ny;
  const size_t nFields = fields.size();
  size_t n_undefined = 0;
  MIUTIL_OPENMP_PARALLEL(fsize, for reduction(+:n_undefined))
  for (int i = 0; i < fsize; i++) {
    int n = 0; // count defined fields
    float m = 0, m2 = 0;
    for (size_t j = 0; j < nFields; j++) {
      if (fieldcalc::is_defined(fDefinedIn[j] == miutil::ALL_DEFINED, fields[j][i], undef)) {
        // see https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Online_algorithm
        const float x = fields[j][i], delta = x - m;
        n += 1;
        m += delta / n;
        m2 += delta * (x - m);
      }
    }
    if (n > 0) {
      fres[i] = sqrt(m2 / n);
    } else {
      fres[i] = undef;
      n_undefined += 1;
    }
  }
  fDefinedOut = miutil::checkDefined(n_undefined, fsize);
  return true;
}

bool extremeValue(int compute, int nx, int ny, const std::vector<float*>& fields, float* fres, ValuesDefined& fDefined, float undef)
{
  // returns min value in each grid point

  //  compute=1 : max - value
  //  compute=2 : min - value
  //  compute=3 : max - index
  //  compute=4 : min - index

  const size_t nFields = fields.size();
  if (nFields == 0)
    return false;
  const int fsize = nx * ny;
  const bool inAllDefined = fDefined == miutil::ALL_DEFINED;
  size_t n_undefined = 0;
  if (compute == 1 || compute == 2) { // max/min value
    MIUTIL_OPENMP_PARALLEL(fsize, for reduction(+:n_undefined))
    for (int i = 0; i < fsize; i++) {
      fres[i] = undef;
      for (size_t j = 0; j < nFields; j++) {
        if (fres[i] == undef || (fieldcalc::is_defined(inAllDefined, fields[j][i], undef) &&
                                 ((compute == 1 && fres[i] < fields[j][i]) || (compute == 2 && fres[i] > fields[j][i])))) {
          fres[i] = fields[j][i];
        }
      }
      if (fres[i] == undef)
        n_undefined += 1;
    }
  } else if (compute == 3 || compute == 4) { // max/min index
    MIUTIL_OPENMP_PARALLEL(fsize, for reduction(+:n_undefined))
    for (int i = 0; i < fsize; i++) {
      fres[i] = undef;
      float tmp = undef;
      for (size_t j = 0; j < nFields; j++) {
        if (tmp == undef ||
            (fieldcalc::is_defined(inAllDefined, fields[j][i], undef) && ((compute == 3 && tmp < fields[j][i]) || (compute == 4 && tmp > fields[j][i])))) {
          tmp = fields[j][i];
          fres[i] = j;
        }
      }
      if (fres[i] == undef)
        n_undefined += 1;
    }
  }
  fDefined = miutil::checkDefined(n_undefined, fsize);
  return true;
}

bool probability(int compute, int nx, int ny, const std::vector<float*>& fields, const std::vector<ValuesDefined>& fDefinedIn, const std::vector<float>& limits,
                 float* fres, ValuesDefined& fDefinedOut, float undef)
{
  // returns probability or number of  value above/below limit in each grid point
  //  compute=1 : above - probability
  //  compute=2 : below - probability
  //  compute=3 : between - probability
  //  compute=4 : above - number
  //  compute=5 : below - number
  //  compute=6 : between - number

  const size_t fsize = nx * ny;
  const size_t lsize = limits.size(), nFields = fields.size();

  const bool check_between = (lsize >= 2) && (compute == 3 || compute == 6);
  const bool check_above = (lsize >= 1) && (compute == 1 || compute == 4 || check_between);
  const bool check_below = (lsize >= 1) && (compute == 2 || compute == 5 || check_between);
  const float value_above = limits[0];
  const float value_below = check_between ? limits[1] : limits[0];

  if (!(check_above || check_below)) {
    MIUTIL_OPENMP_PARALLEL(fsize, for)
    for (size_t i = 0; i < fsize; i++)
      fres[i] = undef;
    fDefinedOut = miutil::NONE_DEFINED;
    return false;
  }

  size_t n_undefined = 0;
  MIUTIL_OPENMP_PARALLEL(fsize, for reduction(+:n_undefined))
  for (size_t i = 0; i < fsize; i++) {
    fres[i] = 0;
    int nfields_defined = 0;
    for (size_t j = 0; j < nFields; j++) {
      if (fDefinedIn[j] != miutil::NONE_DEFINED) {
        nfields_defined += 1;
        const float value = fields[j][i];
        if ((value != undef) && (!check_above || value > value_above) && (!check_below || value < value_below)) {
          fres[i] += 1;
        }
      }
    }

    if (nfields_defined == 0) {
      fres[i] = undef;
      n_undefined += 1;
    } else {
      if (compute < 4)
        fres[i] /= (nfields_defined / 100.0);
    }
  }
  fDefinedOut = miutil::checkDefined(n_undefined, fsize);
  return true;
}

bool neighbourProbFunctions(int nx, int ny, const float* field, const std::vector<float>& constants, int compute, float* fres, ValuesDefined& fDefined,
                            float undef)
{
  // compute=5 : calc probability above (constant[0]=limit)
  // compute=6 : calc probability below (constant[0]=limit)

  if (fDefined != miutil::ALL_DEFINED) {
    return false;
  }

  if (constants.size() < 2) {
    return false;
  }

  int fsize = nx * ny;
  int limit = constants[0];
  int range = constants[1];

  // Calc. probability in each gridcell
  if (compute == 5) {
    MIUTIL_OPENMP_PARALLEL(fsize, for)
    for (int i = 0; i < fsize; i++) {
      fres[i] = field[i] > limit ? 1 : 0;
    }
  } else if (compute == 6) {
    MIUTIL_OPENMP_PARALLEL(fsize, for)
    for (int i = 0; i < fsize; i++) {
      fres[i] = field[i] < limit ? 1 : 0;
    }
  }

  if (range == 0)
    return true;

  // using neighbour method

  // Calc. the summed area table
  float* tmp = new float[nx * ny];
  MIUTIL_OPENMP_PARALLEL(nx, for)
  for (int i = 0; i < nx; i++) {
    tmp[i] = fres[i];
    for (int j = 1; j < ny; j++)
      tmp[i + j * nx] = fres[i + j * nx] + tmp[i + (j - 1) * nx];
  }
  MIUTIL_OPENMP_PARALLEL(ny, for)
  for (int j = 0; j < ny; j++)
    for (int i = 1; i < nx; i++)
      tmp[i + j * nx] += tmp[(i - 1) + j * nx];

  // Calc. the neighbourhood mean for each pixel
  const int N = (2 * range + 1) * (2 * range + 1);
  MIUTIL_OPENMP_PARALLEL(nx-2*range, for)
  for (int i = range; i < nx - range; i++) {
    int imax = i + range;
    for (int j = range; j < ny - range; j++) {
      int jmax = j + range;
      fres[i + j * nx] = tmp[imax + jmax * nx];
      if (i > range) {
        fres[i + j * nx] -= tmp[i - range - 1 + jmax * nx];
        if (j > range)
          fres[i + j * nx] += tmp[i - range - 1 + (j - range - 1) * nx] - tmp[imax + (j - range - 1) * nx];
      } else if (j > range) {
        fres[i + j * nx] -= tmp[imax + (j - range - 1) * nx];
      }
      fres[i + j * nx] /= N;
    }
  }

  delete[] tmp;

  // values can not be calculated close to the border, set values undef
  fDefined = miutil::SOME_DEFINED;
  for (int j = 0; j < range; j++) {
    for (int i = 0; i < nx; i++) {
      fres[i + j * nx] = undef;
    }
  }
  for (int j = range; j < ny - range; j++) {
    for (int i = 0; i < range; i++) {
      fres[i + j * nx] = undef;
    }
    for (int i = nx - range; i < nx; i++) {
      fres[i + j * nx] = undef;
    }
  }
  for (int j = ny - range; j < ny; j++) {
    for (int i = 0; i < nx; i++) {
      fres[i + j * nx] = undef;
    }
  }
  return true;
}

bool neighbourFunctions(int nx, int ny, const float* field, const std::vector<float>& constants, int compute, float* fres, ValuesDefined& fDefined, float undef)
{
  // compute=1 : calc mean value
  // compute=2 : calc max value
  // compute=3 : calc min value
  // compute=4 : calc percentile (constant[0]=percentile)
  // compute=5 : calc probability above (constant[0]=limit)
  // compute=6 : calc probability below (constant[0]=limit)

  if (fDefined != miutil::ALL_DEFINED) {
    return false;
  }

  if (constants.size() < 1 || (constants.size() < 2 && compute > 3)) {
    return false;
  }

  int range = 3, step = 3, limit = 0;
  if (compute < 4) {
    range = constants[0];
    if (constants.size() == 2)
      step = constants[1];
  } else {
    limit = constants[0];
    range = constants[1];
    if (constants.size() == 3)
      step = constants[2];
  }

  if (range > nx || range > ny || range < 1) {
    return false;
  }

  if (step < 1) {
    return false;
  }

  // values can not be calculated close to the border, set values undef
  fDefined = miutil::SOME_DEFINED;
  for (int j = 0; j < range; j++) {
    for (int i = 0; i < nx; i++) {
      fres[i + j * nx] = undef;
    }
  }
  for (int j = range; j < ny - range; j++) {
    for (int i = 0; i < range; i++) {
      fres[i + j * nx] = undef;
    }
    for (int i = nx - range; i < nx; i++) {
      fres[i + j * nx] = undef;
    }
  }
  for (int j = ny - range; j < ny; j++) {
    for (int i = 0; i < nx; i++) {
      fres[i + j * nx] = undef;
    }
  }

  // no. of gridpoints in range
  const float ngridp = miutil::square(2 * range + 1);
  // calc. index corresponding to given percentile and range
  const int ii = ngridp * limit / 100;

  // loop through all gridpoints with given step
  for (int j = range; j < ny - range; j += step) {
    for (int i = range; i < nx - range; i += step) {

      std::vector<float> values;
      float value = 0.0;
      if (compute == 2 || compute == 3) {
        value = field[(i - range) + (j - range) * nx];
      }

      for (int k = j - range; k < j + range + 1; k++) {
        for (int j = i - range; j < i + range + 1; j++) {
          const float thisvalue = field[j + k * nx];
          if (compute == 1) {
            value += thisvalue;
          }
          if ((compute == 2 && thisvalue > value) || (compute == 3 && thisvalue < value)) {
            value = thisvalue;
          } else if (compute == 4) {
            values.push_back(thisvalue);
          }
          if ((compute == 5 && thisvalue > limit) || (compute == 6 && thisvalue < limit)) {
            value++;
          }
        }
      }
      if (compute == 4) {
        std::sort(values.begin(), values.end());
        value = values[ii];
      }
      if (compute == 1 || compute > 4) {
        value /= ngridp;
      }

      // set value in output field
      for (int l = j - (step - 1) / 2; l < j + step / 2 + 1; l++) {
        for (int k = i - (step - 1) / 2; k < i + step / 2 + 1; k++) {
          fres[k + l * nx] = value;
        }
      }
    }
  }
  return true;
}

bool snow_in_cm(int nx, int ny, const float* snow_water, const float* tk2m, const float* td2m, float* snow_cm, ValuesDefined& fDefined, float undef)
{
  /*
Function: snow_in_cm
Parameters: in - snow_water (kg/m2), tk2m (2 m air temperature in kelvin), td2m (2 m dew point temperature in kelvin).
            out - snow depth in cm.
Algorithm:
Computing snow in cm.
 t2m = 2 meter temperature in kelvin
 td = dew point temperature in kelvin

t=(t2m+td)/2

logit(t)=(1-e^(t-274.3)3.5)/(1+e^(t-274.3)3.5 )

mm2cm(t)=0.13/(0.02+0.1((t-252.0)/20)^2 )

fac=logit(t)*mm2cm(t)

snow_cm=smow_water*fac

Reference:
MESAN Mesoskalig analys, SMHI RMK Nr 75, Mars 1997.
   */

  const int fsize = nx * ny;
  bool inAllDefined = fDefined == miutil::ALL_DEFINED;
  size_t n_undefined = 0;
  MIUTIL_OPENMP_PARALLEL(fsize, for reduction(+:n_undefined))
  for (int i = 0; i < fsize; i++) {
    if (fieldcalc::is_defined(inAllDefined, snow_water[i], tk2m[i], td2m[i], undef)) {
      if (snow_water[i] <= 0.) {
        snow_cm[i] = 0.;
        continue;
      }

      float t = (tk2m[i] + td2m[i]) / 2.;

      float logit_t = (1 - exp((t - 274.3) * 3.5)) / (1 + exp((t - 274.3) * 3.5));

      float mm2cm_t = 0.13 / (0.02 + 0.1 * ((t - 252.0) / 20.0) * ((t - 252.0) / 20.0));

      float fac = logit_t * mm2cm_t;
      // We must preserve model consistency
      if (fac <= 1.)
        snow_cm[i] = snow_water[i];
      else
        snow_cm[i] = snow_water[i] * fac;
    } else {
      snow_cm[i] = undef;
      n_undefined += 1;
    }
  }
  fDefined = miutil::checkDefined(n_undefined, fsize);
  return true;
}

} // namespace fieldcalc
} // namespace miutil
