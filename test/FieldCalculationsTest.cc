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


#include <mi_fieldcalc/FieldCalculations.h>
#include <mi_fieldcalc/MetConstants.h>

#include <gtest/gtest.h>

#include <cmath>
#include <iostream>
#include <iomanip>

namespace /* anonymous */ {
const float UNDEF = 12356789, T0 = 273.15;
} // anonymous namespace

namespace {
struct levelhum_params_t {
  int cah; // 'compute' for alevelhum and hlevelhum
  int cp;  // 'compute' for plevelhum
  float t; // temperature input
  float humin; // humidity input
  float p; // pressure input
  float expect; // expected output
  float near; // max deviation from expected
};
} // namespace

TEST(FieldCalculationsTest, absHum)
{
  float out = 2 * UNDEF;
  float t = 293.16;
  float rh = 0.8;
  float exp = 13.82;

  miutil::ValuesDefined fDefined = miutil::ALL_DEFINED;

  EXPECT_TRUE(miutil::fieldcalc::abshum(1, 1, &t, &rh, &out, fDefined, UNDEF));
  EXPECT_NEAR(exp, out, 0.1) << " t=" << t << " rh=" << rh << " abshum=" << out;
  EXPECT_EQ(miutil::ALL_DEFINED, fDefined);
}

TEST(FieldCalculationsTest, XLevelHum)
{
  const levelhum_params_t levelhum_params[] = {
    // alevelhum/hlevelhum and plevelhum have compute numbers >= 5 switched
    {  1,  1, 30.68f+T0, .025f, 1013, 91.9f, 0.1f },
    {  2,  2, 302.71f,   .025f, 1013, 91.9f, 0.1f },
    {  3,  3, 30.68f+T0, 55,    1013, 0.014963f, .000001f },
    {  4,  4, 302.71f,   55,    1013, 0.014963f, .000001f },
    {  5,  7, 30.68f+T0, .015f, 1013, 20.6f, 0.1f },
    {  6,  8, 302.71f,   .015f, 1013, 20.6f, 0.1f },
    {  7,  5, 30.68f+T0, 55,    1013, 20.6f, 0.1f },
    {  8,  6, 302.71f,   55,    1013, 20.6f, 0.1f },
    { -1, -1, 0, 0, 0, 0, 0 }
  };
  const float alevel = 0, blevel = 1;

  for (int i=0; levelhum_params[i].cah >= 0; ++i) {
    const levelhum_params_t& p = levelhum_params[i];

    float humout = 2*UNDEF;
    miutil::ValuesDefined fDefined = miutil::ALL_DEFINED;
    EXPECT_TRUE(miutil::fieldcalc::alevelhum(1, 1, &p.t, &p.humin, &p.p, "celsius", p.cah, &humout, fDefined, UNDEF));
    EXPECT_NEAR(p.expect, humout, p.near) << "alevelhum C i=" << i << " compute=" << p.cah;
    EXPECT_EQ(miutil::ALL_DEFINED, fDefined);

    humout = 2*UNDEF;
    fDefined = miutil::ALL_DEFINED;
    EXPECT_TRUE(miutil::fieldcalc::hlevelhum(1, 1, &p.t, &p.humin, &p.p, alevel, blevel, "celsius", p.cah, &humout, fDefined, UNDEF));
    EXPECT_NEAR(p.expect, humout, p.near) << "hlevelhum C i=" << i << " compute=" << p.cah;
    EXPECT_EQ(miutil::ALL_DEFINED, fDefined);

    humout = 2*UNDEF;
    fDefined = miutil::ALL_DEFINED;
    EXPECT_TRUE(miutil::fieldcalc::plevelhum(1, 1, &p.t, &p.humin, p.p, "celsius", p.cp, &humout, fDefined, UNDEF));
    EXPECT_NEAR(p.expect, humout, p.near) << "plevelhum C i=" << i << " compute=" << p.cp;
    EXPECT_EQ(miutil::ALL_DEFINED, fDefined);

    humout = 2*UNDEF;
    fDefined = miutil::SOME_DEFINED;
    EXPECT_TRUE(miutil::fieldcalc::alevelhum(1, 1, &p.t, &p.humin, &p.p, "celsius", p.cah, &humout, fDefined, UNDEF));
    EXPECT_NEAR(p.expect, humout, p.near) << "alevelhum C i=" << i << " compute=" << p.cah;
    EXPECT_EQ(miutil::ALL_DEFINED, fDefined);

    humout = 2*UNDEF;
    fDefined = miutil::SOME_DEFINED;
    EXPECT_TRUE(miutil::fieldcalc::hlevelhum(1, 1, &p.t, &p.humin, &p.p, alevel, blevel, "celsius", p.cah, &humout, fDefined, UNDEF));
    EXPECT_NEAR(p.expect, humout, p.near) << "hlevelhum C i=" << i << " compute=" << p.cah;
    EXPECT_EQ(miutil::ALL_DEFINED, fDefined);

    humout = 2*UNDEF;
    fDefined = miutil::SOME_DEFINED;
    EXPECT_TRUE(miutil::fieldcalc::plevelhum(1, 1, &p.t, &p.humin, p.p, "celsius", p.cp, &humout, fDefined, UNDEF));
    EXPECT_NEAR(p.expect, humout, p.near) << "plevelhum C i=" << i << " compute=" << p.cp;
    EXPECT_EQ(miutil::ALL_DEFINED, fDefined);

    if (p.cah < 5)
      continue;

    fDefined = miutil::ALL_DEFINED;
    EXPECT_TRUE(miutil::fieldcalc::alevelhum(1, 1, &p.t, &p.humin, &p.p, "kelvin", p.cah, &humout, fDefined, UNDEF));
    EXPECT_NEAR(p.expect+T0, humout, p.near) << "alevelhum K i=" << i << " compute=" << p.cah;
    EXPECT_EQ(miutil::ALL_DEFINED, fDefined);

    fDefined = miutil::ALL_DEFINED;
    EXPECT_TRUE(miutil::fieldcalc::hlevelhum(1, 1, &p.t, &p.humin, &p.p, alevel, blevel, "kelvin", p.cah, &humout, fDefined, UNDEF));
    EXPECT_NEAR(p.expect+T0, humout, p.near) << "hlevelhum K i=" << i << " compute=" << p.cah;
    EXPECT_EQ(miutil::ALL_DEFINED, fDefined);

    fDefined = miutil::ALL_DEFINED;
    EXPECT_TRUE(miutil::fieldcalc::plevelhum(1, 1, &p.t, &p.humin, p.p, "kelvin", p.cp, &humout, fDefined, UNDEF));
    EXPECT_NEAR(p.expect+T0, humout, p.near) << "plevelhum K i=" << i << " compute=" << p.cp;
    EXPECT_EQ(miutil::ALL_DEFINED, fDefined);
  }
}

TEST(FieldCalculationsTest, ALevelTempPerformance)
{
  const float UNDEF = 1e30, T0 = 273.15;
  miutil::ValuesDefined fDefined = miutil::ALL_DEFINED;

  { // compute == 3, T(Kelvin) -> TH
    const int N = 719*929; //1000000;
    const float F = 0.00001;
    float *tk = new float[N], *p = new float[N], *th = new float[N], *expected = new float[N];
    for (int i=0; i<N; ++i) {
      tk[i] = 20 + (i*F) + T0;
      p[i]  = 1005 + (i*F);
      th[i] = 2*UNDEF;
    }
    for (int i=0; i<1; ++i)
      EXPECT_TRUE(miutil::fieldcalc::aleveltemp(1, N, tk, p, "kelvin", 3, th, fDefined, UNDEF));
    for (int i=0; i<N; ++i) {
      const float ex = tk[i] / powf(p[i] * miutil::constants::p0inv, miutil::constants::kappa);
      EXPECT_FLOAT_EQ(ex, th[i]);
    }
    delete[] tk;
    delete[] p;
    delete[] th;
    delete[] expected;
  }
}

namespace {
struct oper_params_t {
  int c; // 'compute'
  float a, b; // input
  miutil::ValuesDefined input_defined;
  float expect; // expected output
};

enum { PLUS = 1, MINUS = 2, MUL = 3, DIV = 4 };

const oper_params_t oper_params[] = {
  { PLUS,  1, 3, miutil::ALL_DEFINED,    4 },
  { PLUS,  1, 3, miutil::SOME_DEFINED, 4,},
  { MINUS, 1, 3, miutil::ALL_DEFINED,    -2 },
  { MINUS, 1, 3, miutil::SOME_DEFINED, -2 },
  { MUL, 1.5, 3, miutil::ALL_DEFINED,    4.5 },
  { MUL, 1.5, 3, miutil::SOME_DEFINED, 4.5 },
  { DIV, 3, 1.5, miutil::ALL_DEFINED,    2 },
  { DIV, 3, 1.5, miutil::SOME_DEFINED, 2 },
  { DIV, 3, 0,   miutil::ALL_DEFINED,    UNDEF },
  { DIV, 3, 0,   miutil::SOME_DEFINED, UNDEF },
  { -1,  0, 0,   miutil::SOME_DEFINED, 0 }
};
} // namespace

TEST(FieldCalculationsTest, XOperX)
{
  float out;
  miutil::ValuesDefined fDefined;

  for (int i=0; oper_params[i].c > 0; ++i) {
    const oper_params_t& op = oper_params[i];

    out = 2*UNDEF;
    fDefined = op.input_defined;
    EXPECT_TRUE(miutil::fieldcalc::fieldOPERfield(op.c, 1, 1, &op.a, &op.b, &out, fDefined, UNDEF));
    EXPECT_EQ(op.expect == UNDEF, miutil::NONE_DEFINED == fDefined) << "fOPf i=" << i;
    EXPECT_NEAR(op.expect, out, 1e-6) << " fOPf i=" << i << " c=" << op.c << " a=" << op.a << " b=" << op.b;

    out = 2*UNDEF;
    fDefined = op.input_defined;
    EXPECT_TRUE(miutil::fieldcalc::fieldOPERconstant(op.c, 1, 1, &op.a, op.b, &out, fDefined, UNDEF));
    EXPECT_EQ(op.expect == UNDEF, miutil::NONE_DEFINED == fDefined) << "fOPc i=" << i;
    EXPECT_NEAR(op.expect, out, 1e-6) << " fOPc i=" << i << " c=" << op.c << " a=" << op.a << " b=" << op.b;

    out = 2*UNDEF;
    fDefined = op.input_defined;
    EXPECT_TRUE(miutil::fieldcalc::constantOPERfield(op.c, 1, 1, op.a, &op.b, &out, fDefined, UNDEF));
    EXPECT_EQ(op.expect == UNDEF, miutil::NONE_DEFINED == fDefined) << "cOPf i=" << i;
    EXPECT_NEAR(op.expect, out, 1e-6) << " cOPf i=" << i << " c=" << op.c << " a=" << op.a << " b=" << op.b;
  }
}

TEST(FieldCalculationsTest, Probability)
{
  const float UNDEF = 123456;

  const int N_ENS = 10;
  float fields[N_ENS];
  std::vector<float*> infields;
  for (size_t i=0; i<N_ENS; ++i)
    infields.push_back(&fields[i]);

  std::fill(fields, fields+N_ENS, UNDEF);
  fields[2] = 940;
  fields[4] = 3500;

  std::vector<miutil::ValuesDefined> defined(N_ENS, miutil::SOME_DEFINED);
  defined[0] = miutil::NONE_DEFINED;
  defined[8] = miutil::NONE_DEFINED;

  std::vector<float> limits(2, 3000);
  float field_out = UNDEF;
  miutil::ValuesDefined defined_out = miutil::NONE_DEFINED;
  EXPECT_TRUE(miutil::fieldcalc::probability(/*compute=*/2 /*probability_below*/,
                                 /*nx=*/1, /*ny=*/1, infields, defined, limits, &field_out, defined_out, UNDEF));
  EXPECT_NEAR(100.0*1/8, field_out, 1e-6);
  EXPECT_EQ(miutil::ALL_DEFINED, defined_out);

  field_out = UNDEF;
  defined_out = miutil::NONE_DEFINED;
  EXPECT_TRUE(miutil::fieldcalc::probability(/*compute=*/1 /*probability_above*/,
                                 /*nx=*/1, /*ny=*/1, infields, defined, limits, &field_out, defined_out, UNDEF));
  EXPECT_NEAR(100.0*1/8, field_out, 1e-6);
  EXPECT_EQ(miutil::ALL_DEFINED, defined_out);

  field_out = UNDEF;
  defined_out = miutil::NONE_DEFINED;
  limits[0] = 4000;
  EXPECT_TRUE(miutil::fieldcalc::probability(/*compute=*/2 /*probability_below*/,
                                 /*nx=*/1, /*ny=*/1, infields, defined, limits, &field_out, defined_out, UNDEF));
  EXPECT_NEAR(100.0*2/8, field_out, 1e-6);
  EXPECT_EQ(miutil::ALL_DEFINED, defined_out);

  field_out = UNDEF;
  defined_out = miutil::NONE_DEFINED;
  limits[0] =  500;
  limits[1] = 4000;
  EXPECT_TRUE(miutil::fieldcalc::probability(/*compute=*/3 /*probability_between*/,
                                 /*nx=*/1, /*ny=*/1, infields, defined, limits, &field_out, defined_out, UNDEF));
  EXPECT_NEAR(100.0*2/8, field_out, 1e-6);
  EXPECT_EQ(miutil::ALL_DEFINED, defined_out);
}

TEST(FieldCalculationsTest, Probability12)
{
  const float UNDEF = 123456;

  const int N_ENS = 10;
  float fields[N_ENS];
  std::vector<float*> infields;
  for (size_t i=0; i<N_ENS; ++i)
    infields.push_back(&fields[i]);

  std::fill(fields, fields+N_ENS, 12);
  fields[3] = fields[5] = UNDEF;

  std::vector<miutil::ValuesDefined> defined(N_ENS, miutil::SOME_DEFINED);

  std::vector<float> limits(1, 3000);
  float field_out = UNDEF;
  miutil::ValuesDefined defined_out = miutil::NONE_DEFINED;
  EXPECT_TRUE(miutil::fieldcalc::probability(/*compute=*/2 /*probability_below*/,
                                 /*nx=*/1, /*ny=*/1, infields, defined, limits, &field_out, defined_out, UNDEF));
  EXPECT_NEAR(80, field_out, 1e-6);
  EXPECT_EQ(miutil::ALL_DEFINED, defined_out);

  field_out = UNDEF;
  defined_out = miutil::NONE_DEFINED;
  EXPECT_TRUE(miutil::fieldcalc::probability(/*compute=*/1 /*probability_above*/,
                                 /*nx=*/1, /*ny=*/1, infields, defined, limits, &field_out, defined_out, UNDEF));
  EXPECT_NEAR(0, field_out, 1e-6);
  EXPECT_EQ(miutil::ALL_DEFINED, defined_out);
}

TEST(FieldCalculationsTest, Neighbour)
{
  const float UNDEF = 123456;

  const int NX = 10, NY = 10;
  float infield[NX * NY], outfield[NX * NY], outfield2[NX * NY];
  std::vector<float> constants;
  miutil::ValuesDefined fDefined = miutil::ALL_DEFINED;
  miutil::ValuesDefined fDefined2 = miutil::ALL_DEFINED;

  constants.push_back(NX+1);
  EXPECT_FALSE(miutil::fieldcalc::neighbourFunctions(NX, NY, infield, constants, 2, outfield, fDefined, UNDEF));

  constants.clear();
  constants.push_back(NX);
  constants.push_back(0);
  EXPECT_FALSE(miutil::fieldcalc::neighbourFunctions(NX, NY, infield, constants, 2, outfield, fDefined, UNDEF));

  std::fill(infield, infield + NX * NY, 0);
  infield[16] = 6;
  std::fill(outfield, outfield + NX * NY, UNDEF / 2);
  constants.clear();
  constants.push_back(3);
  constants.push_back(3);
  EXPECT_TRUE(miutil::fieldcalc::neighbourFunctions(NX, NY, infield, constants, 2, outfield, fDefined, UNDEF));
  EXPECT_EQ(miutil::SOME_DEFINED, fDefined);

  for (int i = 0; i < NX; i++) {
    for (int j = 0; j < NY; j++) {
      float expected;
      if (i < 2 || i >= NX - 2 || j < 2 || j >= NY - 2)
        expected = UNDEF;
      else if (j < 5 )
        expected = 6;
      else
        expected = 0;
      EXPECT_EQ(expected, outfield[i + j*NX]) << "i="<<i<<" j="<<j;
    }
  }

  fDefined = miutil::ALL_DEFINED;
  std::fill(infield, infield + NX * NY, 0);
  infield[25] = 6;
  infield[26] = 6;
  infield[35] = 6;
  infield[36] = 6;
  std::fill(outfield, outfield + NX * NY, UNDEF / 2);
  constants.clear();
  constants.push_back(90);
  constants.push_back(2);
  constants.push_back(1);
  EXPECT_TRUE(miutil::fieldcalc::neighbourFunctions(NX, NY, infield, constants, 4, outfield, fDefined, UNDEF));
  EXPECT_EQ(miutil::SOME_DEFINED, fDefined);

  for (int i = 0; i < NX; i++) {
    for (int j = 0; j < NY; j++) {
      float expected;
      if (i < 2 || i >= NX - 2 || j < 2 || j >= NY - 2)
        expected = UNDEF;
      else if (i > 3 && i < 8 && j > 1 && j < 5 )
        expected = 6;
      else
        expected = 0;
      EXPECT_EQ(expected, outfield[i + j*NX]) << "i="<<i<<" j="<<j;
    }
  }

  fDefined = miutil::ALL_DEFINED;
  fDefined2 = miutil::ALL_DEFINED;
  std::fill(infield, infield + NX * NY, 0);
  infield[25] = 6;
  infield[26] = 6;
  infield[35] = 6;
  infield[36] = 6;
  std::fill(outfield, outfield + NX * NY, UNDEF / 2);
  std::fill(outfield2, outfield2 + NX * NY, UNDEF / 2);
  constants.clear();
  constants.push_back(5);
  constants.push_back(2);
  constants.push_back(1);
  EXPECT_TRUE(miutil::fieldcalc::neighbourFunctions(NX, NY, infield, constants, 5, outfield, fDefined, UNDEF));
  EXPECT_EQ(miutil::SOME_DEFINED, fDefined);
  EXPECT_TRUE(miutil::fieldcalc::neighbourProbFunctions(NX, NY, infield, constants, 5, outfield2, fDefined2, UNDEF));
  EXPECT_EQ(miutil::SOME_DEFINED, fDefined2);

  for (int i = 0; i < NX; i++) {
    for (int j = 0; j < NY; j++) {
      EXPECT_FLOAT_EQ(outfield[i + j * NX], outfield2[i + j * NX]) << "i=" << i << " j=" << j;
    }
  }

  for (int i = 0; i < NX; i++) {
    for (int j = 0; j < NY; j++) {
      float expected;
      if (i < 2 || i >= NX - 2 || j < 2 || j >= NY - 2)
        expected = UNDEF;
      else if (i > 3  && j < 5 )
        expected = 0.16;
      else if (i == 3 && j == 5 )
        expected = 0.04;
      else if (i > 2 && j < 6 )
        expected = 0.08;
      else
        expected = 0;
      EXPECT_EQ(expected, outfield[i + j*NX]) << "i="<<i<<" j="<<j;
    }
  }

  fDefined = miutil::ALL_DEFINED;
  fDefined2 = miutil::ALL_DEFINED;
  std::fill(infield, infield + NX * NY, 0);
  infield[25] = 6;
  infield[26] = 6;
  infield[35] = 6;
  infield[36] = 6;
  std::fill(outfield, outfield + NX * NY, UNDEF / 2);
  std::fill(outfield2, outfield2 + NX * NY, UNDEF / 2);
  constants.clear();
  constants.push_back(5);
  constants.push_back(3);
  constants.push_back(1);
  EXPECT_TRUE(miutil::fieldcalc::neighbourFunctions(NX, NY, infield, constants, 6, outfield, fDefined, UNDEF));
  EXPECT_EQ(miutil::SOME_DEFINED, fDefined);
  EXPECT_TRUE(miutil::fieldcalc::neighbourProbFunctions(NX, NY, infield, constants, 6, outfield2, fDefined2, UNDEF));
  EXPECT_EQ(miutil::SOME_DEFINED, fDefined2);

  for (int i = 0; i < NX; i++) {
    for (int j = 0; j < NY; j++) {
      EXPECT_FLOAT_EQ(outfield[i + j * NX], outfield2[i + j * NX]) << "i=" << i << " j=" << j;
    }
  }

  for (int i = 0; i < NX; i++) {
    for (int j = 0; j < NY; j++) {
      float expected;
      if (i < 3 || i >= NX - 3 || j < 3 || j >= NY - 3)
        expected = UNDEF;
      else if (j < 6)
        expected = 45. / 49.;
      else
        expected = 47. / 49.;
      EXPECT_FLOAT_EQ(expected, outfield[i + j * NX]) << "i=" << i << " j=" << j;
    }
  }
}

TEST(FieldCalculationsTest, ReplaceDefined)
{
  const int n = 2;
  const float in1[n] = { 0, 1 };
  float out1[n] = { -1, -1 };

  miutil::ValuesDefined def1 = miutil::SOME_DEFINED;
  miutil::fieldcalc::replaceDefined(n, 1, in1, 5, out1, def1, 0);
  EXPECT_EQ(0, out1[0]);
  EXPECT_EQ(5, out1[1]);
  EXPECT_EQ(miutil::ALL_DEFINED, def1);

  def1 = miutil::ALL_DEFINED;
  miutil::fieldcalc::replaceDefined(n, 1, in1, 7, out1, def1, -1);
  EXPECT_EQ(7, out1[0]);
  EXPECT_EQ(7, out1[1]);
  EXPECT_EQ(miutil::ALL_DEFINED, def1);

  def1 = miutil::NONE_DEFINED;
  miutil::fieldcalc::replaceDefined(n, 1, in1, 7, out1, def1, -1);
  EXPECT_EQ(-1, out1[0]);
  EXPECT_EQ(-1, out1[1]);
  EXPECT_EQ(miutil::NONE_DEFINED, def1);

  def1 = miutil::SOME_DEFINED;
  miutil::fieldcalc::replaceDefined(n, 1, in1, 1, out1, def1, 1); // value == undef
  EXPECT_EQ(1, out1[0]);
  EXPECT_EQ(1, out1[1]);
  EXPECT_EQ(miutil::NONE_DEFINED, def1);
}

TEST(FieldCalculationsTest, ReplaceUndefined)
{
  const int n = 2;
  const float in1[n] = { 0, 1 };
  float out1[n] = { -1, -1 };

  miutil::ValuesDefined def1 = miutil::SOME_DEFINED;
  miutil::fieldcalc::replaceUndefined(n, 1, in1, 5, out1, def1, 0);
  EXPECT_EQ(5, out1[0]);
  EXPECT_EQ(1, out1[1]);
  EXPECT_EQ(miutil::ALL_DEFINED, def1);

  def1 = miutil::ALL_DEFINED;
  miutil::fieldcalc::replaceUndefined(n, 1, in1, 7, out1, def1, -1);
  EXPECT_EQ(0, out1[0]);
  EXPECT_EQ(1, out1[1]);
  EXPECT_EQ(miutil::ALL_DEFINED, def1);

  def1 = miutil::NONE_DEFINED;
  miutil::fieldcalc::replaceUndefined(n, 1, in1, 7, out1, def1, -1);
  EXPECT_EQ(7, out1[0]);
  EXPECT_EQ(7, out1[1]);
  EXPECT_EQ(miutil::ALL_DEFINED, def1);

  def1 = miutil::SOME_DEFINED;
  miutil::fieldcalc::replaceUndefined(n, 1, in1, 1, out1, def1, 1); // value == undef
  EXPECT_EQ(0, out1[0]);
  EXPECT_EQ(1, out1[1]);
  EXPECT_EQ(miutil::SOME_DEFINED, def1);
}
