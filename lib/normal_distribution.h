#ifndef NORMAL_DISTRIBUTION_H
#define NORMAL_DISTRIBUTION_H

#include <cmath>
#include <utility>

/* Utilities for integrating normal PDF
*/

namespace Integration {
    template <typename T>
    double normal_cdf(T x) {
        return std::erfc(-x / std::sqrt(2)) / 2;
    }

    // Integral of the normal PDF over the cube
    const auto pdf_integral = [](const auto& hc) {
        long double rv = 1;
        for (const auto& [a, b] : hc.intervals)
            rv *= normal_cdf(b) - normal_cdf(a);
        return rv;
    };

    // Minimum and maximum of the normal PDF at the cube
    const auto pdf_minmax = [](const auto& hc) {
        long double min = 0, max = 0;
        constexpr long double common = std::pow(2 * pi(), hc.dimensions / -2.0);

        for (const auto& [a, b]: hc.intervals) {
            if ((a + b) / 2 < 0)
                min += a * a;
            else min += b * b;
        }

        min = common * exp(min / -2);

        for (const auto& [a, b] : hc.intervals) {
            if (a >= 0)
                max += a * a;
            else if (b < 0)
                max += b * b;
        }

        max = common * exp(max / -2);
        return std::make_pair(min, max);
    };
}

#endif // NORMAL_DISTRIBUTION_H
