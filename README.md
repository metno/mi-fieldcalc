# mi-fieldcalc -- MET Norway Meteorological Field Calculations Library

This library provides some operations on meteorological fields for C++
and python (only a subset of functions is accessible from python at
present).

The python binding is implemented using pybind11. In the CMake build,
the python binding may be disabled by setting `ENABLE_PYTHON` to `0`.

Unit tests are implemented using google-test.

## Use with CMake

Build and install and use `FIND_PACKAGE(mi-fieldcalc)`, then link to
`mi-fieldcalc`.
