#ifndef LINEAR_H
#define LINEAR_H

#include <algorithm>
#include <array>
#include <cstddef>
#include <functional>
#include <utility>
#include <vector>

namespace Integration {
    // <vector of coefficients (e), number (d)> pair representing linear equations
    using linear_equation = std::pair<std::vector<long double>, long double>;

    // Maximum of <e, x> + d at hc
    auto linear_max(const auto& hc, const linear_equation::first_type& e,
                    const linear_equation::second_type& d) {
        long double rv = d;
        for (std::size_t i = 0; i < hc.dimensions; i++) {
            const auto& [a, b] = hc.intervals[i];
            rv += e[i] * (e[i] >= 0 ? b : a);
        }
        return rv;
    }

    // Minimum of <e, x> + d at hc
    auto linear_min(const auto& hc, const linear_equation::first_type& e,
                    const linear_equation::second_type& d) {
        long double rv = d;
        for (std::size_t i = 0; i < hc.dimensions; i++) {
            const auto& [a, b] = hc.intervals[i];
            rv += e[i] * (e[i] >= 0 ? a : b);
        }
        return rv;
    }

    // Measure of the cube with one linear restriction <e, x> + d <= 0 applied
    auto single_section_measure(const auto& hc, const linear_equation::first_type& e,
                                const linear_equation::second_type& d) {
        long double min = 0, max = 0, en = 0;
        std::array<long double, hc.dimensions> u, v;
        std::array<bool, hc.dimensions> is_active;
        is_active.fill(true);

        // Calculate the minimum and maximum of <e, x> + d at the cube
        // en is ||<e, e>||
        for (std::size_t i = 0; i < hc.dimensions; i++) {
            const auto& [a, b] = hc.intervals[i];
            if (e[i] >= 0) {
                min += (u[i] = a) * e[i];
                max += (v[i] = b) * e[i];
            } else {
                min += (u[i] = b) * e[i];
                max += (v[i] = a) * e[i];
            }
            en += e[i] * e[i];
        }

        std::function<long double(long double, long double, long double, long double)> f;

        f = [&](long double d, long double min, long double max, long double en) -> long double {
            std::size_t i, active = std::count(is_active.begin(), is_active.end(), true);

            if (max + d <= 0) {
                long double prod = 1;
                for (i = 0; i < is_active.size(); i++) {
                    if (is_active[i]) {
                        const auto& [a, b] = hc.intervals[i];
                        prod *= b - a;
                    }
                }
                return prod;

            } else if (min + d >= 0) {
                return 0;

            } else if (active == 1) {
                for (i = 0; i < is_active.size(); i++)
                    if (is_active[i])
                        break;

                const auto& [a, b] = hc.intervals[i];
                if (e[i] == 0) {
                    return (d > 0) ? 0 : b - a;
                } else {
                    long double u = -d/e[i];
                    if (e[i] < 0) {
                        if (u >= b) return 0;
                        else if (u <= a) return b - a;
                        else return b - u;
                    } else {
                        if (u <= a) return 0;
                        else if (u >= b) return b - a;
                        else return u - a;
                    }
                }

            } else {
                long double t = -(max + d) / en, rv = 0,
                            w, corrected_min, corrected_max, corrected_en;
                for (i = 0; i < is_active.size(); i++) {
                    if (is_active[i]) {
                        const auto& [a, b] = hc.intervals[i];
                        corrected_min = min - u[i] * e[i];
                        corrected_max = max - v[i] * e[i];
                        corrected_en = en - e[i] * e[i];
                        w = v[i] + t * e[i];

                        is_active[i] = false;
                        rv += f(d + a * e[i], corrected_min, corrected_max, corrected_en) * (w - a) / active;
                        rv += f(d + b * e[i], corrected_min, corrected_max, corrected_en) * (b - w) / active;
                        is_active[i] = true;
                    }
                }

                return rv;
            }
        };

        return f(d, min, max, en);

    }
}

#endif // LINEAR_H
