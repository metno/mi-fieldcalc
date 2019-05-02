/* -*- c++ -*-
 mi-fieldcalc - Meteorological Field Calculations

 Copyright (C) 2006-2019 met.no

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

#ifndef MI_FIELDCALC_FIELDCALCULATIONS_H
#define MI_FIELDCALC_FIELDCALCULATIONS_H

#include "FieldDefined.h"

#include <cmath>
#include <string>
#include <vector>

namespace miutil {
namespace fieldcalc {

inline bool is_defined(float in, float undef)
{
  return !std::isnan(in) && in != undef;
}

inline bool is_defined(bool allDefined, float in1, float undef)
{
  return allDefined or is_defined(in1, undef);
}

inline bool is_defined(bool allDefined, float in1, float in2, float undef)
{
  return allDefined or (is_defined(in1, undef) and is_defined(in2, undef));
}

inline bool is_defined(bool allDefined, float in1, float in2, float in3, float undef)
{
  return allDefined or (is_defined(in1, undef) and is_defined(in2, undef) and is_defined(in3, undef));
}

inline bool is_defined(bool allDefined, float in1, float in2, float in3, float in4, float undef)
{
  return allDefined or (is_defined(in1, undef) and is_defined(in2, undef) and is_defined(in3, undef) and is_defined(in4, undef));
}

inline bool is_defined(bool allDefined, float in1, float in2, float in3, float in4, float in5, float undef)
{
  return allDefined or (is_defined(in1, undef) and is_defined(in2, undef) and is_defined(in3, undef) and is_defined(in4, undef) and is_defined(in5, undef));
}

inline bool is_defined(bool allDefined, float in1, float in2, float in3, float in4, float in5, float in6, float undef)
{
  return allDefined or (is_defined(in1, undef) and is_defined(in2, undef) and is_defined(in3, undef) and is_defined(in4, undef) and is_defined(in5, undef) and is_defined(in6, undef));
}

inline bool is_defined(bool allDefined, float in1, float in2, float in3, float in4, float in5, float in6, float in7, float undef)
{
  return allDefined or (is_defined(in1, undef) and is_defined(in2, undef) and is_defined(in3, undef) and is_defined(in4, undef) and is_defined(in5, undef) and is_defined(in6, undef) and is_defined(in7, undef));
}

inline bool is_defined(bool allDefined, float in1, float in2, float in3, float in4, float in5, float in6, float in7, float in8, float undef)
{
  return allDefined or (is_defined(in1, undef) and is_defined(in2, undef) and is_defined(in3, undef) and is_defined(in4, undef) and is_defined(in5, undef) and is_defined(in6, undef) and is_defined(in7, undef) and is_defined(in8, undef));
}

inline bool is_defined(bool allDefined, float in1, float in2, float in3, float in4, float in5, float in6, float in7, float in8, float in9, float undef)
{
  return allDefined or
         (is_defined(in1, undef) and is_defined(in2, undef) and is_defined(in3, undef) and is_defined(in4, undef) and is_defined(in5, undef) and is_defined(in6, undef) and is_defined(in7, undef) and is_defined(in8, undef) and is_defined(in9, undef));
}

inline bool is_defined(bool allDefined, float in1, float in2, float in3, float in4, float in5, float in6, float in7, float in8, float in9, float in10,
                       float undef)
{
  return allDefined or (is_defined(in1, undef) and is_defined(in2, undef) and is_defined(in3, undef) and is_defined(in4, undef) and is_defined(in5, undef) and is_defined(in6, undef) and is_defined(in7, undef) and is_defined(in8, undef) and
                        is_defined(in9, undef) and is_defined(in10, undef));
}

void copy_field(float* fout, const float* fin, size_t fsize);

// parameter order:
// - nx, ny
// - const float* input fields
// - other parameters, "compute" last
// - float* output field
// - fDefined, undef

//---------------------------------------------------
// pressure level (PLEVEL) functions
//---------------------------------------------------

bool pleveltemp(int nx, int ny, const float* tinp, float p, const std::string& unit, int compute, float* tout, ValuesDefined& fDefined, float undef);

bool plevelthe(int nx, int ny, const float* t, const float* rh, float p, int compute, float* the, ValuesDefined& fDefined, float undef);

bool plevelhum(int nx, int ny, const float* t, const float* huminp, float p,
               const std::string& unit, int compute, float* humout, ValuesDefined& fDefined, float undef);

bool pleveldz2tmean(int nx, int ny, const float* z1, const float* z2, float p1, float p2, int compute, float* tmean, ValuesDefined& fDefined, float undef);

bool plevelqvector(int nx, int ny, const float* z, const float* t, const float* xmapr, const float* ymapr, const float* fcoriolis,
                   float p, int compute, float* qcomp, ValuesDefined& fDefined, float undef);

bool plevelducting(int nx, int ny, const float* t, const float* h, float p, int compute, float* duct, ValuesDefined& fDefined, float undef);

bool plevelgwind_xcomp(int nx, int ny, const float* z, const float* xmapr, const float* ymapr, const float* fcoriolis, float* ug, ValuesDefined& fDefined,
                       float undef);

bool plevelgwind_ycomp(int nx, int ny, const float* z, const float* xmapr, const float* ymapr, const float* fcoriolis, float* vg, ValuesDefined& fDefined,
                       float undef);

bool plevelgvort(int nx, int ny, const float* z, const float* xmapr, const float* ymapr, const float* fcoriolis, float* gvort, ValuesDefined& fDefined,
                 float undef);

bool kIndex(int nx, int ny, const float* t500, const float* t700, const float* rh700, const float* t850, const float* rh850,
            float p500, float p700, float p850, int compute, float* kfield, ValuesDefined& fDefined, float undef);

bool ductingIndex(int nx, int ny, const float* t850, const float* rh850, float p850, int compute, float* duct, ValuesDefined& fDefined, float undef);

bool showalterIndex(int nx, int ny, const float* t500, const float* t850, const float* rh850, float p500, float p850,
                    int compute, float* sfield, ValuesDefined& fDefined, float undef);

bool boydenIndex(int nx, int ny, const float* t700, const float* z700, const float* z1000, float p700, float p1000,
                 int compute, float* bfield, ValuesDefined& fDefined, float undef);

bool sweatIndex(int nx, int ny, const float* t850, const float* t500, const float* td850, const float* td500, const float* u850, const float* v850,
                const float* u500, const float* v500, float* sindex, ValuesDefined& fDefined, float undef);

//---------------------------------------------------
// hybrid model level (HLEVEL) functions
//---------------------------------------------------

bool hleveltemp(int nx, int ny, const float* tinp, const float* ps, float alevel, float blevel,
                const std::string& unit, int compute, float* tout, ValuesDefined& fDefined, float undef);

bool hlevelthe(int nx, int ny, const float* t, const float* q, const float* ps, float alevel, float blevel, int compute, float* the, ValuesDefined& fDefined,
               float undef);

bool hlevelhum(int nx, int ny, const float* t, const float* huminp, const float* ps, float alevel, float blevel, const std::string& unit, int compute, float* humout,
               ValuesDefined& fDefined, float undef);

bool hlevelducting(int nx, int ny, const float* t, const float* h, const float* ps, float alevel, float blevel, int compute, float* duct,
                   ValuesDefined& fDefined, float undef);

bool hlevelpressure(int nx, int ny, const float* ps, float alevel, float blevel, float* p, ValuesDefined& fDefined, float undef);

//---------------------------------------------------
// atmospheric model level (ALEVEL) functions
//---------------------------------------------------

bool aleveltemp(int nx, int ny, const float* tinp, const float* p, const std::string& unit, int compute, float* tout, ValuesDefined& fDefined, float undef);

bool alevelthe(int nx, int ny, const float* t, const float* q, const float* p, int compute, float* the, ValuesDefined& fDefined, float undef);

bool alevelhum(int nx, int ny, const float* t, const float* huminp, const float* p,
               const std::string& unit, int compute, float* humout, ValuesDefined& fDefined, float undef);

bool alevelducting(int nx, int ny, const float* t, const float* h, const float* p, int compute, float* duct, ValuesDefined& fDefined, float undef);

//---------------------------------------------------
// isentropic level (ILEVEL) function
//---------------------------------------------------

bool ilevelgwind(int nx, int ny, const float* mpot, const float* xmapr, const float* ymapr, const float* fcoriolis, float* ug, float* vg,
                 ValuesDefined& fDefined, float undef);

//---------------------------------------------------
// ocean depth level (OZLEVEL) functions
//---------------------------------------------------

bool seaSoundSpeed(int nx, int ny, const float* t, const float* s, float z, int compute, float* soundspeed, ValuesDefined& fDefined, float undef);

//---------------------------------------------------
// level (pressure) independant functions
//---------------------------------------------------

bool cvtemp(int nx, int ny, const float* tinp, int compute, float* tout, ValuesDefined& fDefined, float undef);

bool cvhum(int nx, int ny, const float* t, const float* huminp, const std::string& unit, int compute, float* humout, ValuesDefined& fDefined, float undef);

bool abshum(int nx, int ny, const float* t, const float* rhum, float* abshumout, ValuesDefined& fDefined, float undef);

bool vectorabs(int nx, int ny, const float* u, const float* v, float* ff, ValuesDefined& fDefined, float undef);

bool relvort(int nx, int ny, const float* u, const float* v, const float* xmapr, const float* ymapr, float* rvort, ValuesDefined& fDefined, float undef);

bool absvort(int nx, int ny, const float* u, const float* v, const float* xmapr, const float* ymapr, const float* fcoriolis, float* avort,
             ValuesDefined& fDefined, float undef);

bool divergence(int nx, int ny, const float* u, const float* v, const float* xmapr, const float* ymapr, float* diverg, ValuesDefined& fDefined, float undef);

bool advection(int nx, int ny, const float* f, const float* u, const float* v, const float* xmapr, const float* ymapr, float hours, float* advec,
               ValuesDefined& fDefined, float undef);

bool gradient(int nx, int ny, const float* field, const float* xmapr, const float* ymapr, int compute, float* fgrad, ValuesDefined& fDefined, float undef);

bool shapiro2_filter(int nx, int ny, float* field, float* fsmooth, ValuesDefined& fDefined, float undef);

bool windCooling(int nx, int ny, const float* t, const float* u, const float* v, int compute, float* dtcool, ValuesDefined& fDefined, float undef);

bool underCooledRain(int nx, int ny, const float* precip, const float* snow, const float* tk, float precipMin, float snowRateMax,
                     float tcMax, float* undercooled, ValuesDefined& fDefined, float undef);

bool thermalFrontParameter(int nx, int ny, const float* t, const float* xmapr, const float* ymapr, float* tfp, ValuesDefined& fDefined, float undef);

bool pressure2FlightLevel(int nx, int ny, const float* pressure, float* flightlevel, ValuesDefined& fDefined, float undef);

bool momentumXcoordinate(int nx, int ny, const float* v, const float* xmapr, const float* fcoriolis, float fcoriolisMin, float* mxy, ValuesDefined& fDefined,
                         float undef);

bool momentumYcoordinate(int nx, int ny, const float* u, const float* ymapr, const float* fcoriolis, float fcoriolisMin, float* nxy, ValuesDefined& fDefined,
                         float undef);

bool jacobian(int nx, int ny, const float* field1, const float* field2, const float* xmapr, const float* ymapr, float* fjacobian, ValuesDefined& fDefined,
              float undef);

bool vesselIcingOverland(int nx, int ny, const float* airtemp, const float* seatemp, const float* u, const float* v, const float* sal, const float* aice,
                         float* icing, ValuesDefined& fDefined, float undef);

bool vesselIcingMertins(int nx, int ny, const float* airtemp, const float* seatemp, const float* u, const float* v, const float* sal, const float* aice,
                        float* icing, ValuesDefined& fDefined, float undef);

bool vesselIcingModStall(int nx, int ny, const float* sal, const float* wave, const float* x_wind, const float* y_wind, const float* airtemp, const float* rh,
                         const float* sst, const float* p, const float* Pw, const float* aice, const float* depth, const float vs,
                         const float alpha, const float zmin, const float zmax, float* icing, ValuesDefined& fDefined, float undef);

bool vesselIcingMincog(int nx, int ny, const float* sal, const float* wave, const float* x_wind, const float* y_wind, const float* airtemp, const float* rh,
                       const float* sst, const float* p, const float* Pw, const float* aice, const float* depth, const float vs,
                       const float alpha, const float zmin, const float zmax, const int alt, float* icing, ValuesDefined& fDefined, float undef);

bool values2classes(int nx, int ny, const float* fvalue, float* fclass, const std::vector<float>& values, ValuesDefined& fDefined, float undef);

void minvalueFields(int nx, int ny, const float* field1, const float* field2, float* fres, ValuesDefined& fDefined, float undef);

void minvalueFieldConst(int nx, int ny, const float* field1, const float value, float* fres, ValuesDefined& fDefined, float undef);

void maxvalueFields(int nx, int ny, const float* field1, const float* field2, float* fres, ValuesDefined& fDefined, float undef);

void maxvalueFieldConst(int nx, int ny, const float* field1, const float value, float* fres, ValuesDefined& fDefined, float undef);

void absvalueField(int nx, int ny, const float* field, float* fres, ValuesDefined& fDefined, float undef);

void log10Field(int nx, int ny, const float* field, float* fres, ValuesDefined& fDefined, float undef);

void pow10Field(int nx, int ny, const float* field, float* fres, ValuesDefined& fDefined, float undef);

void logField(int nx, int ny, const float* field, float* fres, ValuesDefined& fDefined, float undef);

void expField(int nx, int ny, const float* field, float* fres, ValuesDefined& fDefined, float undef);

void powerField(int nx, int ny, const float* field, float value, float* fres, ValuesDefined& fDefined, float undef);

void replaceUndefined(int nx, int ny, const float* field, float value, float* fres, ValuesDefined& fDefined, float undef);

void replaceDefined(int nx, int ny, const float* field, float value, float* fres, ValuesDefined& fDefined, float undef);

bool fieldOPERfield(int compute, int nx, int ny, const float* field1, const float* field2, float* fres, ValuesDefined& fDefined, float undef);

bool fieldOPERconstant(int compute, int nx, int ny, const float* field, float value, float* fres, ValuesDefined& fDefined, float undef);

bool constantOPERfield(int compute, int nx, int ny, float value, const float* field, float* fres, ValuesDefined& fDefined, float undef);

bool sumFields(int nx, int ny, const std::vector<float*>& fields, float* fres, ValuesDefined& fDefined, float undef);

bool meanValue(int nx, int ny, const std::vector<float*>& fields, const std::vector<ValuesDefined>& fDefinedIn, float* fres, ValuesDefined& fDefinedOut,
               float undef);

bool stddevValue(int nx, int ny, const std::vector<float*>& fields, const std::vector<ValuesDefined>& fDefinedIn, float* fres, ValuesDefined& fDefinedOut,
                 float undef);

bool extremeValue(int compute, int nx, int ny, const std::vector<float*>& fields, float* fres, ValuesDefined& fDefined, float undef);

bool probability(int compute, int nx, int ny, const std::vector<float*>& fields, const std::vector<ValuesDefined>& fDefinedIn, const std::vector<float>& limits,
                 float* fres, ValuesDefined& fDefinedOut, float undef);

bool neighbourProbFunctions(int nx, int ny, const float* field, const std::vector<float>& constants, int compute, float* fres, ValuesDefined& fDefined,
                            float undef);

bool neighbourFunctions(int nx, int ny, const float* field, const std::vector<float>& constants, int compute, float* fres, ValuesDefined& fDefined,
                        float undef);

bool snow_in_cm(int nx, int ny, const float* snow_water, const float* tk2m, const float* td2m, float* snow_cm, ValuesDefined& fDefined, float undef);

} // namespace fieldcalc
} // namespace miutil

#endif // MI_FIELDCALC_FIELDCALCULATIONS_H
