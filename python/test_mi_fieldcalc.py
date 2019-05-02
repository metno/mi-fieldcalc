# Metlibs-FieldCalc, python/test_mifieldcalc.py
#
# Copyright (C) 2019 met.no
#
# Contact information:
# Norwegian Meteorological Institute
# Box 43 Blindern
# 0313 OSLO
# NORWAY
# email: diana@met.no
#
# This library is free software; you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation; either version 2.1 of the License, or
# (at your option) any later version.
#
# This library is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
# or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
# License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
# USA.

import mi_fieldcalc as fc
import numpy as np
import os.path
import unittest

test_srcdir = os.path.join(os.path.dirname(__file__), "..", "..", "test")

class TestFieldCalc(unittest.TestCase):

    def test_abshum(self):
        tk = np.array([[293.16]])
        rh = np.array([[0.8]])
        ah = fc.abshum(tk, rh, -1)
        self.assertIsNotNone(ah)
        self.assertAlmostEqual(13.83, ah[0,0], delta=0.02)

if __name__ == '__main__':
    unittest.main()
