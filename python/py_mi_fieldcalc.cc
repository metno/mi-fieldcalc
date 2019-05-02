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

py::object py_abshum(npa_t t, npa_t rhum)
{
  if (t.ndim() != 2)
    return py::none();

  miutil::ValuesDefined fDefined = miutil::SOME_DEFINED;
  const float undef = miutil::UNDEF;
  const int nx = t.shape(0), ny = t.shape(1);
  npa_t abshum({nx, ny});

  if (!(same_dims(t, rhum) && same_dims(t, abshum)))
    return py::none();

  if (!miutil::fieldcalc::abshum(nx, ny, t.data(), rhum.data(), abshum.mutable_data(), fDefined, undef))
    return py::none();

  return std::move(abshum);
}

py::object py_vesselIcingOverland(npa_t airtemp, npa_t seatemp, npa_t u, npa_t v, npa_t sal, npa_t aice, float undef)
{
#define SAME_DIMS(array) same_dims(airtemp, array)
  if (!(SAME_DIMS(seatemp) && SAME_DIMS(u) && SAME_DIMS(v) && SAME_DIMS(sal) && SAME_DIMS(aice)))
    return py::none();
#undef SAME_DIMS

  miutil::ValuesDefined fDefined = miutil::SOME_DEFINED;
  const int nx = airtemp.shape(0), ny = airtemp.shape(1);
  const py::array::ShapeContainer shp(airtemp.shape(), airtemp.shape() + airtemp.ndim());
  py::array_t<float> icing(shp);
  if (miutil::fieldcalc::vesselIcingOverland(nx, ny, airtemp.data(), seatemp.data(), u.data(), v.data(), sal.data(), aice.data(),
                                             icing.mutable_data(), fDefined, undef))
    return std::move(icing);
  else
    return py::none();
}

py::object py_vesselIcingMertins(npa_t airtemp, npa_t seatemp, npa_t u, npa_t v, npa_t sal, npa_t aice, float undef)
{
#define SAME_DIMS(array) same_dims(airtemp, array)
  if (!(SAME_DIMS(seatemp) && SAME_DIMS(u) && SAME_DIMS(v) && SAME_DIMS(sal) && SAME_DIMS(aice)))
    return py::none();
#undef SAME_DIMS

  miutil::ValuesDefined fDefined = miutil::SOME_DEFINED;
  const int nx = airtemp.shape(0), ny = airtemp.shape(1);
  const py::array::ShapeContainer shp(airtemp.shape(), airtemp.shape() + airtemp.ndim());
  py::array_t<float> icing(shp);
  if (miutil::fieldcalc::vesselIcingMertins(nx, ny, airtemp.data(), seatemp.data(), u.data(), v.data(), sal.data(), aice.data(),
                                            icing.mutable_data(), fDefined, undef))
    return std::move(icing);
  else
    return py::none();
}

py::object py_vesselIcingModStall(npa_t sal, npa_t wave, npa_t x_wind, npa_t y_wind, npa_t airtemp, npa_t rh, npa_t sst, npa_t p, npa_t Pw, npa_t aice,
                                  npa_t depth, float vs, float alpha, float zmin, float zmax, float undef)
{
#define SAME_DIMS(array) same_dims(sal, array)
  if (!(SAME_DIMS(wave) && SAME_DIMS(x_wind) && SAME_DIMS(y_wind) && SAME_DIMS(airtemp) && SAME_DIMS(rh) && SAME_DIMS(sst) && SAME_DIMS(p) && SAME_DIMS(Pw) &&
        SAME_DIMS(aice) && SAME_DIMS(depth)))
    return py::none();
#undef SAME_DIMS

  miutil::ValuesDefined fDefined = miutil::SOME_DEFINED;
  const int nx = sal.shape(0), ny = sal.shape(1);
  const py::array::ShapeContainer shp(sal.shape(), sal.shape() + sal.ndim());
  py::array_t<float> icing(shp);
  if (miutil::fieldcalc::vesselIcingModStall(nx, ny, sal.data(), wave.data(), x_wind.data(), y_wind.data(), airtemp.data(), rh.data(), sst.data(), p.data(), Pw.data(),
                                             aice.data(), depth.data(), icing.mutable_data(), vs, alpha, zmin, zmax, fDefined, undef))
    return std::move(icing);
  else
    return py::none();
}

py::object py_vesselIcingMincog(npa_t sal, npa_t wave, npa_t x_wind, npa_t y_wind, npa_t airtemp, npa_t rh, npa_t sst, npa_t p, npa_t Pw, npa_t aice,
                                npa_t depth, float vs, float alpha, float zmin, float zmax, float alt, float undef)
{
#define SAME_DIMS(array) same_dims(sal, array)
  if (!(SAME_DIMS(wave) && SAME_DIMS(x_wind) && SAME_DIMS(y_wind) && SAME_DIMS(airtemp) && SAME_DIMS(rh) && SAME_DIMS(sst) && SAME_DIMS(p) && SAME_DIMS(Pw) &&
        SAME_DIMS(aice) && SAME_DIMS(depth)))
    return py::none();
#undef SAME_DIMS

  miutil::ValuesDefined fDefined = miutil::SOME_DEFINED;
  const int nx = sal.shape(0), ny = sal.shape(1);
  const py::array::ShapeContainer shp(sal.shape(), sal.shape() + sal.ndim());
  py::array_t<float> icing(shp);
  if (miutil::fieldcalc::vesselIcingMincog(nx, ny, sal.data(), wave.data(), x_wind.data(), y_wind.data(), airtemp.data(), rh.data(), sst.data(), p.data(), Pw.data(),
                                           aice.data(), depth.data(), icing.mutable_data(), vs, alpha, zmin, zmax, alt, fDefined, undef))
    return std::move(icing);
  else
    return py::none();
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

  m.def("abshum", py_abshum);

  m.def("vesselIcingOverland", py_vesselIcingOverland);
  m.def("vesselIcingMertins", py_vesselIcingMertins);
  m.def("vesselIcingModStall", py_vesselIcingModStall);
  m.def("vesselIcingMincog", py_vesselIcingMincog);
}
