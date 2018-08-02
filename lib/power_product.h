#ifndef POWER_PRODUCT_H
#define POWER_PRODUCT_H

#include <algorithm>
#include <vector>
#include <utility>
#include <cstddef>
#include <cmath>

namespace Integration {
    // Represents a function of the form x1^a1 * x2^a2 ... xN^aN,
    // i.e. a monomial/power product
    struct power_product {
        std::vector<unsigned int> exponents;

        power_product(std::vector<unsigned int> _exponents)
                      : exponents(std::move(_exponents)) { }

        // Integral of the function over the cube
        template <typename Cube>
        long double integral(const Cube& hc) const {
            long double rv = 1;
            for (std::size_t i = 0; i < hc.dimensions; i++) {
                const auto& [a, b] = hc.intervals[i];
                const unsigned int exponent = exponents[i] + 1;
                rv *= (std::pow(b, exponent) - std::pow(a, exponent)) / exponent;
            }
            return rv;
        }

        // Minimum and maximum of the function at the cube
        template <typename Cube>
        std::pair<long double, long double> minmax(const Cube& hc) const {
            std::vector<long double> values{1};
            long double a, b, first, second;

            for (std::size_t i = 0; i < exponents.size(); i++) {
                std::vector<long double> temp;
                temp.reserve(values.size() * 3);

                std::tie(a, b) = hc.intervals[i];
                // NOTE: Assumption that 0 ^ p = 0 for any p
                first = (a == 0) ? 0 : std::pow(a, exponents[i]);
                second = (b == 0) ? 0 : std::pow(b, exponents[i]);

                for (const auto& v : values) {
                    temp.push_back(v * first);
                    temp.push_back(v * second);
                    if (exponents[i] % 2 == 0 && a <= 0 && b >= 0) {
                        // 0 ^ exponents[i]
                        temp.push_back(0);
                    }
                }

                values = std::move(temp);
            }

            const auto& [min, max] = std::minmax_element(values.begin(), values.end());
            return {*min, *max};

        }
    };
}

#endif // POWER_PRODUCT_H
