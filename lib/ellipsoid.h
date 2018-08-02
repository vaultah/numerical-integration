#ifndef ELLIPSOID_H
#define ELLIPSOID_H

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <memory>
#include <vector>
#include <queue>
#include <utility>
#include <tuple>
#include "const.h"
#include "integrationresult.h"
#include "linear.h"

/* Utilities for integrating over regions restricted by an equation
 * of the form  a1(x1 - c1)^2 + a2(x2 - c2)^2  + ... + aN(xN - cN)^2 - d <= 0
*/

namespace Integration {
    class ellipsoid
    {
    private:
        std::vector<long double> coeffs, center;
        long double d;
        // Sum of coefficients
        long double sum;

    public:
        ellipsoid(std::vector<long double> _coeffs, std::vector<long double> _center, long double _d)
                  : coeffs(std::move(_coeffs)), center(std::move(_center)), d(_d) {
            sum = 0;
            for (const auto& x : coeffs)
                sum += x;
        }

        template <typename Cube>
        std::pair<long double, long double> measure_estimates(const Cube& hc) const {
            long double gc = -d, tau = 0, a, b, medium;
            std::vector<long double> e(hc.dimensions);

            for (std::size_t i = 0; i < hc.dimensions; i++) {
                std::tie(a, b) = hc.intervals[i];
                medium = (a + b) / 2;
                e[i] = 2 * coeffs[i] * (medium - center[i]);
                gc += coeffs[i] * std::pow(medium - center[i], 2) - e[i] * medium;
                tau = std::max(tau, b - a);
            }

            return {single_section_measure(hc, e, gc + sum * tau * tau / 4),
                    single_section_measure(hc, e, gc)};
        }

        template <typename Cube>
        REGION_STATE contains(const Cube& hc) const {
            long double temp, a, b;
            std::size_t i;

            temp = 0;
            for (i = 0; i < hc.dimensions; i++) {
                std::tie(a, b) = hc.intervals[i];
                temp += coeffs[i] * std::pow((a + b < 2 * center[i] ? a : b) - center[i], 2);
            }

            if (temp - d <= 0)
                return CONTAINED;

            temp = 0;
            for (i = 0; i < hc.dimensions; i++) {
                std::tie(a, b) = hc.intervals[i];
                temp += coeffs[i] * std::pow(std::clamp(center[i], a, b) - center[i], 2);
            }

            if (temp - d >= 0)
                return REJECTED;

            return INDEFINITE;
        }

        template <typename Cube, typename Integral, typename Function>
        Result<Cube> integrate(const Cube& cube, Integral Int, Function f,
                               unsigned max_splits, bool return_cubes = false) const {
            REGION_STATE state;
            auto result = Result<Cube>{};
            auto hc = std::make_shared<Cube>(cube);
            result.origin = hc;
            unsigned int depth = 0;

            using cube_info = std::pair<std::shared_ptr<Cube>, unsigned int>;
            std::queue<cube_info> cubes {{ { hc, depth } }};

            for ( ; !cubes.empty(); cubes.pop()) {
                std::tie(hc, depth) = cubes.front();
                state = contains(*hc);

                switch (state) {
                    case INDEFINITE:
                        if (depth >= max_splits) {
                            const auto& [mlow, mhigh] = measure_estimates(*hc);
                            const auto& [flow, fhigh] = f(*hc);
                            result.sum += std::min(0.0L,  flow) * mhigh + std::max(0.0L,  flow) * mlow;
                            result.error += ((std::max(0.0L, fhigh) - std::min(0.0L,  flow)) * mhigh +
                                             (std::min(0.0L, fhigh) - std::max(0.0L,  flow)) * mlow);
                        } else {
                            depth++;
                            for (auto p : hc->split())
                                cubes.push({ p, depth });
                            continue;
                        }
                        break;

                    case CONTAINED:
                        result.sum += Int(*hc);
                        break;

                    default:
                        break;
                }

                // Destroy the cube unless return_cubes is true
                if (return_cubes)
                    result.cubes.push_back({ hc, state });
            }

            return result;
        }
    };
}

#endif // ELLIPSOID_H
