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
inline constexpr T pow2(T x)
{
  return x * x;
}

template <typename T>
inline T icing_f1(T t)
{
  return T(0.6112) * std::exp(T(17.67) * t / (t + T(243.5)));
}

template <typename T>
inline T icing_f10MK(T t, T M, T K)
{
  return (M - T(0.2) * t) - K * 10 * icing_f1(t);
}

template <typename T>
inline T kT4(T T_celsius)
{
  const T sigma = 5.67e-8;               // Stefan-Boltzmanns const.
  return sigma * pow2(pow2(T_celsius + miutil::constants::t0));
}

} // namespace

namespace miutil {
namespace fieldcalc {

bool vesselIcingOverland(int nx, int ny, const float* airtemp, const float* seatemp, const float* u, const float* v, const float* sal, const float* aice,
                         float* icing, ValuesDefined& fDefined, float undef)
{
  // Based on: Overland (1990)
  // All temperatures in degrees celsius

  const int fsize = nx * ny;

  const double A = 2.73e-2;
  const double B = 2.91e-4;
  const double C = 1.84e-6;

  const bool inAllDefined = fDefined == miutil::ALL_DEFINED;
  size_t n_undefined = 0;
  MIUTIL_OPENMP_PARALLEL(fsize, for reduction(+:n_undefined))
  for (int i = 0; i < fsize; i++) {
    if (fieldcalc::is_defined(inAllDefined, airtemp[i], seatemp[i], u[i], v[i], sal[i], aice[i], undef) && aice[i] < 0.4) {
      /* Freezing point of sea water from Stallabrass (1980) in Celcius*/
      const double Tf = (-0.002 - 0.0524 * sal[i]) - 6.0E-5 * pow2(sal[i]);

      if (seatemp[i] < Tf) {
        icing[i] = undef;
        n_undefined = false;
      } else {
        const double ff = miutil::absval(u[i], v[i]);
        const double ppr = ff * (Tf - airtemp[i]) / (1 + 0.3 * (seatemp[i] - Tf));
        icing[i] = A * ppr + B * (ppr * ppr) + C * ppr * ppr * ppr;
      }
    } else {
      icing[i] = undef;
      n_undefined += 1;
    }
  }
  fDefined = miutil::checkDefined(n_undefined, fsize);
  return true;
}

bool vesselIcingMertins(int nx, int ny, const float* airtemp, const float* seatemp, const float* u, const float* v, const float* sal, const float* aice,
                        float* icing, ValuesDefined& fDefined, float undef)
{
  // Based on: H.O. Mertins : Icing on fishing vessels due to spray, Marine Observer No.221, 1968
  // All temperatures in degrees celsius

  const int fsize = nx * ny;
  const bool inAllDefined = fDefined == miutil::ALL_DEFINED;
  size_t n_undefined = 0;
  MIUTIL_OPENMP_PARALLEL(fsize, for reduction(+:n_undefined))
  for (int i = 0; i < fsize; i++) {
    if (fieldcalc::is_defined(inAllDefined, airtemp[i], seatemp[i], u[i], v[i], sal[i], aice[i], undef) && aice[i] < 0.4) {
      /* Freezing point of sea water from Stallabrass (1980) in Celcius*/
      const double Tf = (-0.002 - 0.0524 * sal[i]) - 6.0E-5 * (sal[i] * sal[i]);

      if (seatemp[i] < Tf) {
        icing[i] = undef;
        n_undefined += 1;
      } else {
        const double ff = miutil::absval(u[i], v[i]);
        const double temperature = airtemp[i];
        const double sst = seatemp[i];
        if (ff >= 10.8) {
          double temp1, temp2, temp3;
          if (ff < 17.2) {
            temp1 = -1.15 * sst - 4.3;
            temp2 = -1.5 * sst - 10;
            temp3 = -10000;
          } else if (ff < 20.8) {
            temp1 = -0.6 * sst - 3.2;
            temp2 = -1.05 * sst - 5.6;
            temp3 = -1.75 * sst - 12.5;
          } else if (ff < 28.5) {
            temp1 = -0.3 * sst - 2.6;
            temp2 = -0.66 * sst - 3.32;
            temp3 = -1.325 * sst - 7.651;
          } else {
            temp1 = -0.14 * sst - 2.28;
            temp2 = -0.3 * sst - 2.6;
            temp3 = -1.16 * sst - 5.22;
          }

          if (temperature > -2) {
            icing[i] = 0;
          } else if (temperature > temp1) {
            icing[i] = 0.8333;
          } else if (temperature > temp2) {
            icing[i] = 2.0833;
          } else {
            if (temperature <= temp3 || ff < 17.2) {
              icing[i] = 4.375;
            } else {
              icing[i] = 6.25;
            }
          }
        } else {
          icing[i] = 0;
        }
      }
    } else {
      icing[i] = undef;
      n_undefined += 1;
    }
  }
  fDefined = miutil::checkDefined(n_undefined, fsize);
  return true;
}

bool vesselIcingModStall(int nx, int ny, const float* sal, const float* wave, const float* x_wind, const float* y_wind, const float* airtemp, const float* rh,
                         const float* sst, const float* p, const float* Pw, const float* aice, const float* depth, float* icing, const float vs,
                         const float alpha, const float zmin, const float zmax, ValuesDefined& fDefined, float undef)
{

  // Modified Stallabrass (described in Henry (1995), Samuelsen et.al. (2015))
  // All temperatures in degrees celsius

  const int fsize = nx * ny;

  const double num = zmax - zmin;
  const int number = (num * 2 + 1);

  if (zmax < zmin || fmod(num, 1) != 0) {
    return false;
  }

  if (vs < 0 || alpha < 0 || zmin < 0 || zmax < 0) {
    return false;
  }

  const bool inAllDefined = fDefined == miutil::ALL_DEFINED;
  size_t n_undefined = 0;
  MIUTIL_OPENMP_PARALLEL(fsize, for reduction(+:n_undefined))
  for (int i = 0; i < fsize; i++) {

    if (fieldcalc::is_defined(inAllDefined, sal[i], wave[i], x_wind[i], y_wind[i], airtemp[i], rh[i], sst[i], p[i], aice[i], depth[i], undef) &&
        aice[i] < 0.4) {
      /* Program to calculate freezing seaspray. Modified Stallabrass.  */
      /* From Brown (1991,2011) */
      /* Equation to be solved. */
      /* Ri = (Rw*cw*(Ts-Td))/lf + ha*(Ts-Ta)/lf + he/lf * (esat(Ts) - rh*esat(Ta)); */
      /* n=Ri/Rw; */
      /* Ts = (1+n)*Tf; n<1. Temperature of freezing brine. */

      /*Calculate VR from vs and angle */
      double c = (9.81 / (2 * M_PI)) * Pw[i];

      // Calculating c in shallow water.
      if (depth[i] <= c * Pw[i] && c != 0) {
        // Setting start error and start velocity
        c = 1.0;
        double err = 1.0;
        double c_new = 0;
        int j = 0;
        while (err > 1e-5) {
          c_new = (9.81 * Pw[i] / (2 * M_PI)) * tanh(2 * M_PI * depth[i] / (Pw[i] * c));
          err = std::abs(c_new - c);
          c = c_new;
          j = j + 1;
          if (j > 10000) {
            c = 0.0;
            break;
          }
        }
      }

      double Vr = c - vs * cos(alpha);

      /* Calculate v from x_wind and y_wind */
      double v = miutil::absval(x_wind[i], y_wind[i]);

      /* Freezing point of sea water from Stallabrass (1980) */
      double Tf = (-0.002 - 0.0524 * sal[i]) - 6.0E-5 * (sal[i] * sal[i]);

      /*  eq. 2.15 Henry, 1995, based on 3m cylinder -5 deg. */
      const double ha_ = 5.17, ha = ha_ * pow(v, 0.8);

      /*  eq. 2.16 Henry, 1995, based on 3m cylinder -5 deg.  */
      double ratio = 89.5 / ha_; // originally: * pow(v, 0.8) / ha

      /* Find water droplet temperature from eq: */
      /*  dTd/dt = 0.2 * (t-Td) *(1+0.622*(lv/P*CP)*(ea-etd)/(t-Td)); */

      double tau = 11.25 - v / 4.0;

      double k1 = sst[i];

      /* Spray residence time from Zakrewski (1986) p.44 */
      /* Low tau gives Td=sst */
      if (tau > 0.0) {
        double K = 311000.0 / ((p[i] / 10.0) * 1005.0);
        double M = 0.2 * airtemp[i] + K * rh[i] * icing_f1(airtemp[i]);
        double h = tau / 50.0;

        double y = sst[i];

        /* Use Runge Kutta method */
        for (int counter = 0; counter < 50; counter++) {

          k1 = (M - 0.2 * y) - K * icing_f1(y);
          double y2 = y + 0.5 * h * k1;
          double k2 = (M - 0.2 * y2) - K * icing_f1(y2);
          double y3 = y + 0.5 * h * k2;
          y2 = (M - 0.2 * y3) - K * icing_f1(y3);
          double y4 = y + h * y2;
          y += h * ((1.0 / 6.0) * (((k1 + 2.0 * k2) + 2.0 * y2) + ((M - 0.2 * y4) - K * icing_f1(y4))));
          k1 = y;
        }
      }

      /* section A.1.2 Stallabrass, 1980 */
      /* section B2 Stallabrass, 1980 */
      /* Iterative method */
      /* From Ross Brown (1991,2011) */
      double ice = 0;
      for (int counter = 0; counter < number; counter++) {
        /* Liquid water content, from Zakrweski spray cloud (1987) */
        /*  Spray flux: eq. 2.2 Henry, 1995 */
        double rw = 6.46E-5 * wave[i] * (Vr * Vr) * exp(-0.55 * (zmin + 0.5 * counter)) * v;

        /* Setting N=0 in expression as a start */
        double N = 0.0;
        double N1 = N;

        /* Setting starterror */
        double err = 1.0;

        /* Running loop when error>=tolerance */
        int j = 0;
        double Ts = Tf;
        double ri = 0;
        while (err >= 1.0E-5 && N >= 0 && N <= 1) {
          Ts = (1.0 + N) * Tf;
          ri = (0.012012012 * rw * (Ts - k1) + (ha / 333000.0) * ((Ts - airtemp[i]) + ratio * (icing_f1(Ts) - rh[i] * icing_f1(airtemp[i]))));
          N1 = (ri / rw);
          err = fabs(N1 - N);
          N = N1;
          j = j + 1;
          if (j > 1000) {
            N = 0.0;
            break;
          }
        }

        /* set n=0 for negative values */
        if (N < 0.0) {
          N = 0.0;
        } else if (N > 1.0) { /* Dry icing mode */
          N = 1.0;
        }
        /* cm/hr */
        ice += N * (rw / 890.0) * 3600.0 * 100.0;

      } /*end calculation over all z*/

      icing[i] = std::abs(ice / number);

    } else {
      icing[i] = undef;
      n_undefined += 1;
    }
  }
  fDefined = miutil::checkDefined(n_undefined, fsize);
  return true;
}

template <typename T>
struct FreezeFracZero {
  T Sw, Ta, ha, he, ea, RH, rw, Tsp, Lwdown, Swdown;
  T operator()(T N) const;
};

template <typename T>
T FreezeFracZero<T>::operator()(T N) const
{
  T cw = 4000;
  T lfs = 3.33e5 * 0.7;
  T Sb = Sw / (1 - N * (1 - 0.3));
  T Ts = T(-54.1126) * (Sb / (1000 - Sb));
  T es = 10 * icing_f1(Ts);
  T Qc = ha * (Ts - Ta);
  T Qe = he * (es - RH * ea);
  T Qd = rw * cw * (Ts - Tsp);
  T Lwup = kT4(Ts);
  T Qr = Lwup - Lwdown - 0.44 * Swdown;
  T ri = (1 / lfs) * (Qc + Qe + Qd + Qr);
  T N1 = (ri / rw);
  return N1 - N;
}

template <typename T>
struct IcingF10MK {
  T M, K;
  T operator()(T t) const { return icing_f10MK(t, M, K); }
};

template <typename T>
inline int sgn(T d)
{
  return d < 0 ? -1 : (d > 0 ? +1 : 0); // -1, 0, +1. FIXED
}

template<typename V>
inline bool same_sign(V a, V b)
{
  return (a>0) == (b>0);
}

template<typename Func, typename V>
V bisection(const Func& func, V a, V b, V EPSILON)
{
  V ffa = func(a);
  V ffb = func(b);
  if ((ffa > 0) == (ffb > 0)) {
    return 0;
  }

  const int iterations_max = 100;
  const int iterations = std::min(static_cast<int>(std::log2((b - a) / EPSILON)), iterations_max);
  V c;
  int j = 0;
  for (; j < iterations; ++j) {
    // Find middle point
    c = (a + b) / 2;

    // Check if middle point is root
    const V ffc = func(c);
    if (ffc == 0)
      return c;

    // Decide the side to repeat the steps
    if (!same_sign(ffc, ffa)) {
      b = c;
      ffb = ffc;
    } else {
      a = c;
      ffa = ffc;
    }
  }
  if (j >= iterations_max)
    c = 0;
  return c;
}

template<typename Func, typename V>
V regula_falsi(const Func& func, V a, V b, V EPSILON)
{
  const int iterations_max = 100;
  int side=0;

  V ffa = func(a);
  V ffb = func(b);
  V c = a;
  for (int j = 0; j < iterations_max; j++) {
    c = (ffa*b - ffb*a) / (ffa - ffb);
    if (std::abs(b-a) < EPSILON * std::abs(b+a))
      break;
    V ffc = func(c);
    if (same_sign(ffc, ffb)) {
      b = c;
      ffb = ffc;
      if (side==-1)
        ffa /= 2;
      side = -1;
    } else if (same_sign(ffa, ffc)) {
      a = c;
      ffa = ffc;
      if (side==+1)
        ffb /= 2;
      side = +1;
    } else {
      break;
    }
  }
  return c;
}

template<typename Func, typename V>
V runge_kutta(const Func& func, V y, V tau)
{
  const int steps = 50;
  const V h = tau / steps, h2 = h/2;
  for (int s = 0; s < steps; s++) {
    V k1 = h2*func(y);
    V k2 = h *func(y + k1);
    V k3 = h *func(y + k2/2);
    V k4 = h2*func(y + k3);
    y += (k1 + k2 + k3 + k4)/3;
  }
  return y;
}

template<typename V>
V vesselIcingMincog(const V sal, const V wave, const V x_wind, const V y_wind, const V airtemp, const V rh,
                        const V sst, const V p, const V Pw, /*const V aice,*/ const V depth, const V vs,
                        const V alpha, const V zmin, const V zmax, const int alt)
{
  // MINCOG described in Samuelsen et.al. (2017) with two options.
  // alt == 1: MINCOG org
  // alt == 2: MINCOG adj. (ref. PhD thesis Samuelsen (2017))
  // If not used as input parameter:
  // const V alt=1; //1, MINCOG, Other: MINCOG adj.
  // All temperatures in degrees celsius

  /* Calculate v from x_wind and y_wind */
  const V v = miutil::absval(x_wind, y_wind);
  if (v < 1 || wave < 0.1) {
    // no icing if v<1 m/s or wave <0.1 m
    return 0;
  }

  /* Program to calculate freezing seaspray. MINCOG (Marine Icing model for the Norwegian COast Guard).  */
  // Samuelsen et al. (2017)
  // Equation to be solved.
  // Qf = Qc + Qe + Qd + Qr

  /*Calculate VR from vs and angle */
  const V c_0 = 9.81 / (2 * M_PI) * Pw; // Wave speed
  V c = c_0;

  // Calculating c in shallow water.
  if (depth <= c * Pw && c_0 != 0) {
    // Setting start error and start velocity
    c = 1;
    int j = 0;
    const V a = 2 * M_PI * depth / Pw;
    for (; j<1000; ++j) {
      const V c_new = c_0 * tanh(a / c);
      const V err = std::abs(c_new - c);
      c = c_new;
      if (err <= 1e-5)
        break;
    }
    if (j >= 1000)
      c = 0;
  }

  const V cos_alpha = cos(alpha);
  const V Vr = c - vs * cos_alpha; // Relative speed between waves and boat

  // Period between spraying. Wave length divided by relative wave speed
  const V tper = std::abs(c * Pw / Vr);
  if (tper <= 0)
    return 0;

  const V beta = alpha, sin_beta = sin(beta);

  // 13 March 2018. Changed Wr calculation.
  const V Wrx = std::abs(v * cos(beta) - vs);
  const V Wry = std::abs(v * sin_beta);
  const V Wr_inv = 1/miutil::absval(Wrx, Wry);

  // 13 March 2018. Changed ha:
  const V hax = 6.0617 * std::pow(Wrx, 1.82);
  const V hay = 4.8496 * std::pow(Wry, 1.8);
  const V ha = (hax + hay) / (Wrx + Wry);

  // 13 March 2018.
  const V tdur = 0.1230 + 0.7008 * std::abs(Vr * wave) / std::max<V>(v, 5);

  // Spray frequency. New expr. from 13 March 2018.
  const V Nf = 1 / (4 * tper); // Spray freq. for coast guard vessel

  /* Find water droplet temperature from eq: */
  /*  dTd/dt = 0.2 * (t-Td) *(1+0.622*(lv/P*CP)*(ea-etd)/(t-Td)); */

  // New calc. from 13 March 2018.
  // Simplified trajectory calculation
  // Find distance according to relative angle
  const V beta_r = M_PI - std::asin(v * sin_beta * Wr_inv);
  V cos_beta_r, cos_2_beta_r, sin_beta_r_2;
  if (beta_r <= (M_PI / 2)) {
    const V br = 91 * M_PI / 180;
    sin_beta_r_2 = pow2(std::sin(br));
    cos_beta_r = std::cos(br);
    cos_2_beta_r = std::cos(2*br);
  } else if (beta_r > (M_PI)) {
    const V br = M_PI;
    sin_beta_r_2 = pow2(std::sin(br));
    cos_beta_r = std::cos(br);
    cos_2_beta_r = std::cos(2*br);
  } else {
    const V br = beta_r;
    sin_beta_r_2 = pow2(std::sin(br));
    cos_beta_r = std::cos(br);
    cos_2_beta_r = std::cos(2*br);
  }

  // Large ellipse fitting the curve of the perimeter of the front part of KV Nordkapp seen from above. Samuelsen et al. (2017)
  const V r0 = 13.18;
  const V a0 = 32.88;
  const V b0 = 6.605;
  const V a0_2 = pow2(a0), b0_2 = pow2(b0), r0_2 = pow2(r0);
  const V c0 = std::sqrt(2) * a0 * b0 * std::sqrt((b0_2 - a0_2) * cos_2_beta_r + a0_2 + b0_2 - 2 * r0_2 * sin_beta_r_2);
  const V r = (r0 * 2 * b0_2 * cos_beta_r + c0) / ((b0_2 - a0_2) * cos_2_beta_r + a0_2 + b0_2);
  // Assuming that the droplets follow a straight line from initial position in coordinate system not following boat. In reality the droplets follows a
  // curved trajectory.

  // Without drag fource. Calculations in MINCOG shows that drag force slows time approximately to the V of straight line trajectory without drag tau
  // = (r/Wr)
  const V tau_const = r * Wr_inv;
  const V beta_deg = beta * (180 / M_PI);
  const V drag = -0.0046 * beta_deg + 2.1912; // Simplified interpolation between traj. model and model without drag! Mean value between v = 1, 4, 15 and 30 m/s.
  const V tau = tau_const * drag;   // drag = t_flight/tau_const.

  const V ea = 10 * icing_f1(airtemp);

  const V K = 0.2 * 0.622 * 2.5E6 / (p * 1005.0);
  const V M = 0.2 * airtemp + K * rh * ea;
  const IcingF10MK<V> icingF10MK = { M, K };
  const V y = runge_kutta(icingF10MK, (V)sst, tau);

  // 13 March 2018. Spray-cloud temperature mean value between individual droplet temperature and SST.
  const V Td = y;
  const V Tsp = 0.5 * (Td + sst);

  // Component of droplet velocity hitting plate.
  const V Vdz = 6.67;                           // Terminal velocity for droplet dr = 0.002 m.
  const V Vdcomp = Wrx * 0.9962 + Vdz * 0.0872; // Plate on KV Nordkapp tilting 5 degrees from the vertical

  V lwc0;
  if (alt == 1) {
    lwc0 = 6.36E-5 * wave * pow2(Vr); // MINCOG org.
  } else {
    // MINCOG adj method.
    const V lambda = c * Pw, dl = 4*M_PI*depth/lambda;
    const V cg = (c / 2) * (1 + dl / std::sinh(dl));
    const V Vgr = cg - vs * cos_alpha;
    lwc0 = 9.5205E-4 * pow2(wave) * std::sqrt(wave / lambda) * Vgr;
  }
  lwc0 = std::abs(lwc0);

  // Setting start salinity
  const V he = ha * 1738.6 / p;
  const V Ta = airtemp;
  const V Vf = (1 + std::cos(85 * M_PI / 180)) / 2; // View factor

  // Simplification not using radiation from model.
  const V Swdown_model = 0; // May be commented out.
  const V eps_atm = 0.7;    // Emissivity of atmosphere.
  // V Lwdown = 0;
  const V Lwdown = eps_atm * kT4(airtemp); // May be replaced by model
  const V Swdown = Swdown_model * Vf;

  float icing = 0;
  const V num = zmax - zmin;
  const int number = (num * 2 + 1);
  for (int counter = 0; counter < number; counter++) {
    /* Liquid water content, from Zakrweski spray cloud (1987), Samuelsen et al. (2015 */
    const V lwc = lwc0 * std::exp(-0.55 * (zmin + 0.5 * counter));

    /*  Spray flux: eq. 2.2 Henry (1995), Samuelsen et al. (2017) */
    const V rw = lwc * Vdcomp * Nf * tdur;

    const FreezeFracZero<V> ffz = { sal, Ta, ha, he, ea, rh, rw, Tsp, Lwdown, Swdown };

    V initial_guess1 = -0.5;
    V initial_guess2 = 1.3;

#if 0
    // Aitken accelarator
    const V x1 = 0; // N=0
    const V x2 = ffz(x1) + x1;
    const V x3 = ffz(x2) + x2;
    const V delx1 = x2 - x1;
    const V delx2 = x3 - x2;
    const V R = x3 - pow2(delx2) / (delx2 - delx1);

    V X = R, Y = R;
    int mysign0 = sgn(ffz(X));
    int mysign1 = mysign0; // Y == X sgn(ffz(Y));
    const V dx = 0.1;

    // Finding interval where zero point exists.
    for(int j1 = 0; j1 < 1000 && X <= 7.5; ++j1) {
      X += dx;
      const int nextsign0 = sgn(ffz(X));
      if (nextsign0 != mysign0) {
        initial_guess1 = X - dx;
        initial_guess2 = X;
        break;
      }
      Y -= dx;
      const int nextsign1 = sgn(ffz(Y));
      if (nextsign1 != mysign1) {
        initial_guess1 = Y + dx;
        initial_guess2 = Y;
        break;
      }
      mysign0 = nextsign0;
      mysign1 = nextsign1;
    }

    const V N = regula_falsi(ffz, initial_guess1, initial_guess2, V(1e-5));
#else
    const V N = bisection(ffz, initial_guess1, initial_guess2, V(1e-5));
#endif

    icing += rw * constrain_value(N, V(0), V(1));
  } // end calculation over all z

  return std::abs(icing / number) * V(3600.0 * 100.0 / 890.0); // cm/hr
}

bool vesselIcingMincog(int nx, int ny, const float* sal, const float* wave, const float* x_wind, const float* y_wind, const float* airtemp, const float* rh,
                       const float* sst, const float* p, const float* Pw, const float* aice, const float* depth, float* icing, const float vs,
                       const float alpha, const float zmin, const float zmax, const int alt, ValuesDefined& fDefined, float undef)
{
  // MINCOG described in Samuelsen et.al. (2017) with two options.
  // alt == 1: MINCOG org
  // alt == 2: MINCOG adj. (ref. PhD thesis Samuelsen (2017))
  // If not used as input parameter:
  // const double alt=1; //1, MINCOG, Other: MINCOG adj.
  // All temperatures in degrees celsius

  if (vs < 0 || alpha < 0 || zmin < 0 || zmax < 0 || zmax < zmin || fmod(zmax - zmin, 1) != 0)
    return false;

  const int fsize = nx * ny;
  const bool inAllDefined = fDefined == miutil::ALL_DEFINED;
  size_t n_undefined = 0;
  MIUTIL_OPENMP_PARALLEL(fsize, for reduction(+:n_undefined))
  for (int i = 0; i < fsize; i++) {
    if (fieldcalc::is_defined(inAllDefined, sal[i], wave[i], x_wind[i], y_wind[i], airtemp[i], rh[i], sst[i], p[i], aice[i], depth[i], undef) && aice[i] < 0.4) {
      icing[i] = vesselIcingMincog(sal[i], wave[i], x_wind[i], y_wind[i], airtemp[i], rh[i], sst[i], p[i], Pw[i], /*aice[i],*/ depth[i], vs, alpha, zmin, zmax, alt);
    } else {
      icing[i] = undef;
      n_undefined += 1;
    }
  }
  fDefined = miutil::checkDefined(n_undefined, fsize);
  return true;
}

} // namespace fieldcalc
} // namespace miutil
