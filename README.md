ClausenFunctions.jl
==========

[![test](https://github.com/Expander/ClausenFunctions.jl/actions/workflows/build.yml/badge.svg)](https://github.com/Expander/ClausenFunctions.jl/actions/workflows/build.yml)

The ClausenFunctions.jl package provides Julia implementations of the
Clausen functions.


Example
-------

```.jl
using ClausenFunctions

cl1(1.0)     # Standard Clausen function Cl_1(x)
cl2(1.0)     # Standard Clausen function Cl_2(x)
cl3(1.0)     # Standard Clausen function Cl_3(x)
cl4(1.0)     # Standard Clausen function Cl_4(x)
cl5(1.0)     # Standard Clausen function Cl_5(x)
cl6(1.0)     # Standard Clausen function Cl_6(x)
cl(10, 1.0)  # Standard Clausen function Cl_n(x) for integer n > 0
sl(10, 1.0)  # Glaisher-Clausen function Sl_n(x) for integer n > 0
```


Documentation
-------------

[https://docs.juliahub.com/ClausenFunctions/](https://docs.juliahub.com/ClausenFunctions/)


Copying
-------

ClausenFunctions.jl is licenced under the MIT License.
