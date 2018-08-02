#include <iostream>
#include <functional>
#include "hypercube.h"
#include "normal_distribution.h"
#include "polygon.h"
#include "2d_viewer.h"

using namespace Integration;
using namespace std;
using namespace placeholders;

int main() {
    HyperCube<2> cube({{ {-5, 5}, {-5, 5} }});
    polygon poly({ {{1, 1}, -4}, {{-3, 1}, -5}, {{1, -2}, -6} });
    auto splits = 6;
    auto result = poly.integrate(cube, pdf_integral, pdf_minmax, splits, true);
    cout << "Low estimate: " << result.sum << std::endl;
    cout << "High estimate: " << result.sum + result.error << std::endl;
    cout << "Error: " << result.error << std::endl;
    view_result(result);
    return 0;
}
