# Numerical integration

This repo contains implementations of some of numerical integration methods described by Nefedov V.N. in his paper
*"Regarding approximation of multidimensional integrals with specified precision"*, 1991
[\[cut and altered copy, PDF, Russian\]](https://drive.google.com/file/d/1RlCIVFFe37wZgHAGqaxhaZq0HWm_6dCr/view).

That includes facilities to integrate PDF of multivariate normal distribution and monomials over polygons and ellipses
without restrictions on the number of variables. Although, for performance reasons, the algorithms and implementation become
impractical in cases where the number of variables is large (≥ 5).

The code in the repository requires a compiler with C++17 support and Qt libraries installed
(for viewing integration region approximations, which only works for two-dimensional cases anyway).

# Example

Given these input parameters:

 - Bounding box: *-5 ≤ x ≤ 5*, *-5 ≤ y ≤ 5*
 - Integration region: { *x + y - 4 ≤ 0*, *-3x + y - 5 ≤ 0*, *x - 2y - 6 ≤ 0* }
 - Function: PDF of bivariate normal distribution

you could use the following code to calculate the low estimate for the integral and the error
(controlled by the *maximum allowed number of cube partitions* -- here `max_splits`)

```c++
using namespace Integration;

HyperCube<2> cube(-5, 5);
polygon poly({ {{1, 1}, -4}, {{-3, 1}, -5}, {{1, -2}, -6} });
unsigned int max_splits = ...;

auto result = poly.integrate(cube, pdf_integral, pdf_minmax, max_splits);
std::cout << "Low estimate: " << result.sum << std::endl;
std::cout << "Error: " << result.error << std::endl;
```

Here's the table of results for values of `max_splits` from 2 to 7:

| `max_splits` | Low estimate | Error       |
| ------------ | ------------ | ----------- |
| 2            | 0.00607456   | 3.26703     |
| 3            | 0.79278      | 0.428885    |
| 4            | 0.897573     | 0.085253    |
| 5            | 0.926847     | 0.019397    |
| 6            | 0.934527     | 0.00450675  |
| 7            | 0.93645      | 0.0010915   |

and below is the GIF of images of region approximations for these cases.
The interior of the region consists of light blue squares, while the boundary of the region consists
of dark blue squares

![](https://i.imgur.com/IKsz23W.gif)

The basic idea is that as `max_splits` increases, the boundary of the region becomes more accurately
approximated, and thus the error decreases about 4 times with each increment of `max_splits` by one.
