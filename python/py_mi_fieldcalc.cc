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

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <mi_fieldcalc/FieldCalculations.h>

namespace py = pybind11;
namespace fc = miutil::fieldcalc;

namespace {

typedef py::array_t<float, py::array::c_style | py::array::forcecast> npa_t;

bool same_dims(const npa_t& a1, const npa_t& a2)
{
  const ssize_t a1_ndim = a1.ndim();
  if (a1_ndim != a2.ndim())
    return false;
  const ssize_t* a1_shp = a1.shape();
  const ssize_t* a2_shp = a2.shape();
  for (ssize_t d = 0; d < a1_ndim; ++d)
    if (a1_shp[d] != a2_shp[d])
      return false;
  return true;
}

template<typename T>
inline bool same_dims(const npa_t&, T)
{
  return true;
}

template <typename T2, typename... TN>
bool same_dims(npa_t var1, T2 var2, TN... var3)
{
  return same_dims(var1, var2) && same_dims(var1, var3 ...);
}

inline const float* extract(npa_t npa) { return npa.data(); }

template<typename T>
inline T extract(const T& t) { return t; }

template <typename F, typename T1, typename ... TN>
py::object py_wrap_2d(F f, float undef, T1 a1, TN ... an)
{
  if (a1.ndim() != 2 || !same_dims(a1, an ...))
    return py::none();

  const py::array::ShapeContainer shp1(a1.shape(), a1.shape() + a1.ndim());
  py::array_t<float> result(shp1);

  const int nx = a1.shape(0), ny = a1.shape(1);
  miutil::ValuesDefined fDefined = miutil::SOME_DEFINED;

  bool ok;
  {
    py::gil_scoped_release release;
    ok = f(nx, ny, extract(a1), extract(an) ..., result.mutable_data(), fDefined, undef);
  }
  if (!ok)
    return py::none();

  return std::move(result);
}

py::object py_kIndex(npa_t t500, npa_t t700, npa_t rh700, npa_t t850, npa_t rh850,
            float p500, float p700, float p850, int compute, float undef)
{
  return py_wrap_2d(fc::kIndex, undef, t500,  t700,  rh700,  t850,  rh850, p500, p700, p850, compute);
}

py::object py_ductingIndex(npa_t t850, npa_t rh850, float p850, int compute, float undef)
{
  return py_wrap_2d(fc::ductingIndex, undef, t850, rh850, p850, compute);
}

py::object py_showalterIndex(npa_t t500, npa_t t850, npa_t rh850, float p500, float p850, int compute, float undef)
{
  return py_wrap_2d(fc::showalterIndex, undef, t500, t850, rh850, p500, p850, compute);
}

py::object py_boydenIndex(npa_t t700, npa_t z700, npa_t z1000, float p700, float p1000, int compute, float undef)
{
  return py_wrap_2d(fc::boydenIndex, undef, t700, z700, z1000, p700, p1000, compute);
}

py::object py_sweatIndex(npa_t t850, npa_t t500, npa_t td850, npa_t td500, npa_t u850, npa_t v850,
                npa_t u500, npa_t v500, float undef)
{
  return py_wrap_2d(fc::sweatIndex, undef, t850, t500, td850, td500, u850, v850, u500, v500);
}

py::object py_seaSoundSpeed(npa_t t, npa_t s, float z, int compute, float undef)
{
  return py_wrap_2d(fc::seaSoundSpeed, undef, t, s, z, compute);
}

py::object py_cvtemp(npa_t tinp, int compute, float undef)
{
  return py_wrap_2d(fc::cvtemp, undef, tinp, compute);
}

py::object py_cvhum(npa_t t, npa_t huminp, const std::string& unit, int compute, float undef)
{
    return py_wrap_2d(fc::cvhum, undef, t, huminp, unit, compute);
}

py::object py_abshum(npa_t t, npa_t rhum, float undef)
{
  return py_wrap_2d(fc::abshum, undef, t, rhum);
}

py::object py_windCooling(npa_t t, npa_t u, npa_t v, int compute, float undef)
{
  return py_wrap_2d(fc::windCooling, undef, t, u, v, compute);
}

py::object py_underCooledRain(npa_t precip, npa_t snow, npa_t tk, float precipMin, float snowRateMax, float tcMax, float undef)
{
  return py_wrap_2d(fc::underCooledRain, undef, precip, snow, tk, precipMin, snowRateMax, tcMax);
}

py::object py_vesselIcingOverland(npa_t airtemp, npa_t seatemp, npa_t u, npa_t v, npa_t sal, npa_t aice, float undef)
{
  return py_wrap_2d(fc::vesselIcingOverland, undef, airtemp, seatemp, u, v, sal, aice);
}

py::object py_vesselIcingMertins(npa_t airtemp, npa_t seatemp, npa_t u, npa_t v, npa_t sal, npa_t aice, float undef)
{
  return py_wrap_2d(fc::vesselIcingMertins, undef, airtemp, seatemp, u, v, sal, aice);
}

py::object py_vesselIcingModStall(npa_t sal, npa_t wave, npa_t x_wind, npa_t y_wind, npa_t airtemp, npa_t rh, npa_t sst, npa_t p, npa_t Pw, npa_t aice,
                                  npa_t depth, float vs, float alpha, float zmin, float zmax, float undef)
{
  return py_wrap_2d(fc::vesselIcingModStall, undef, sal, wave, x_wind, y_wind, airtemp, rh, sst, p, Pw, aice, depth, vs, alpha, zmin, zmax);
}

py::object py_vesselIcingMincog(npa_t sal, npa_t wave, npa_t x_wind, npa_t y_wind, npa_t airtemp, npa_t rh, npa_t sst, npa_t p, npa_t Pw, npa_t aice,
                                npa_t depth, float vs, float alpha, float zmin, float zmax, int alt, float undef)
{
  return py_wrap_2d(fc::vesselIcingMincog, undef, sal, wave, x_wind, y_wind, airtemp, rh, sst, p, Pw, aice, depth, vs, alpha, zmin, zmax, alt);
}

} // namespace

PYBIND11_MODULE(mi_fieldcalc, m)
{
  // clang-format off
  py::enum_<miutil::ValuesDefined>(m, "ValuesDefined")
          .value("ALL_DEFINED", miutil::ALL_DEFINED)
          .value("NONE_DEFINED", miutil::NONE_DEFINED)
          .value("SOME_DEFINED", miutil::SOME_DEFINED);
  // clang-format on

#define M_DEF_FC(x) m.def(#x, py_ ## x)
  M_DEF_FC(kIndex);
  M_DEF_FC(ductingIndex);
  M_DEF_FC(showalterIndex);
  M_DEF_FC(boydenIndex);
  M_DEF_FC(sweatIndex);

  M_DEF_FC(seaSoundSpeed);

  M_DEF_FC(cvtemp);
  M_DEF_FC(cvhum);
  M_DEF_FC(abshum);

  M_DEF_FC(windCooling);
  M_DEF_FC(underCooledRain);

  M_DEF_FC(vesselIcingOverland);
  M_DEF_FC(vesselIcingMertins);
  M_DEF_FC(vesselIcingModStall);
  M_DEF_FC(vesselIcingMincog);
}
