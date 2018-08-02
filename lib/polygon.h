#ifndef POLYGON_H
#define POLYGON_H

#include <cstddef>
#include <algorithm>
#include <vector>
#include <functional>
#include <numeric>
#include <utility>
#include <queue>
#include <tuple>
#include "const.h"
#include "linear.h"

/* Utilities for integrating over regions restricted by
 * equations of the form <e, x> + d <= 0 (linear), i.e. polygons
 */

namespace Integration {
    class polygon
    {
    public:
        std::vector<linear_equation> equations;

        polygon(std::vector<linear_equation> _equations)
                : equations(std::move(_equations)) { }

        template <typename Cube, typename Integral, typename Function>
        auto integrate(const Cube& cube, Integral Int, Function f,
                       unsigned max_splits, bool return_cubes = false) const {
            REGION_STATE state;
            Result<Cube> result{};
            unsigned int depth = 0;
            auto hc = std::make_shared<Cube>(cube);
            result.origin = hc;

            std::vector<std::size_t> all_equations(equations.size());
            std::iota(all_equations.begin(), all_equations.end(), 0);

            using cube_info = std::tuple<std::shared_ptr<Cube>, unsigned int, std::vector<std::size_t>>;
            std::queue<cube_info> cubes {{ {hc, depth, std::move(all_equations)} }};

            for ( ; !cubes.empty(); cubes.pop()) {
                auto& triplet = cubes.front();
                hc = std::get<0>(triplet);
                depth = std::get<1>(triplet);
                std::vector<std::size_t> boundaries = std::move(std::get<2>(triplet));

                boundaries.erase(
                    std::remove_if(boundaries.begin(), boundaries.end(),
                        [&](const auto& i) {
                            const auto& [e, d] = equations[i];
                            return linear_max(*hc, e, d) <= 0;
                        }
                    ), boundaries.end());

                if (boundaries.empty()) {
                    state = CONTAINED;
                } else {
                    state = INDEFINITE;

                    for (const auto& i : boundaries) {
                        const auto& [e, d] = equations[i];
                        if (linear_min(*hc, e, d) >= 0) {
                            state = REJECTED;
                            break;
                        }
                    }
                }

                switch (state) {
                    case INDEFINITE:
                        if (depth >= max_splits) {
                            long double mlow, mhigh, flow, fhigh;
                            if (boundaries.size() > 1) {
                                mlow = 0;
                                mhigh = hc->volume();
                            } else {
                                const auto& [e, d] = equations[boundaries[0]];
                                mlow = mhigh = single_section_measure(*hc, e, d);
                            }

                            std::tie(flow, fhigh) = f(*hc);
                            result.sum += std::min(0.0L,  flow) * mhigh + std::max(0.0L,  flow) * mlow;
                            result.error += ((std::max(0.0L, fhigh) - std::min(0.0L,  flow)) * mhigh +
                                             (std::min(0.0L, fhigh) - std::max(0.0L,  flow)) * mlow);
                        } else {
                            depth++;
                            for (auto p : hc->split())
                                cubes.push({ p, depth, boundaries });
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

#endif // POLYGON_H
