#ifndef HYPERCUBE_H
#define HYPERCUBE_H

#include <array>
#include <algorithm>
#include <utility>
#include <vector>
#include <cstddef>
#include <cassert>
#include <iostream>
#include <memory>
#include <typeinfo>
#include "const.h"
#include "integrationresult.h"


namespace Integration {
    template<unsigned int N, typename T = double>
    class HyperCube {
    public:
        static const unsigned int dimensions = N;
        using interval_type = std::pair<T, T>;
        using intervals_type = std::array<interval_type, N>;
        using parts_type = std::array<std::shared_ptr<HyperCube>, 1 << N>;

        intervals_type intervals;

        HyperCube(T a, T b) {
            intervals.fill(interval_type{a, b});
        }

        // Defaults to unit cube
        HyperCube() : HyperCube(0, 1) { }

        HyperCube(intervals_type pairs) : intervals(std::move(pairs)) { }

        // Split the current cube into 2^N smaller cubes
        auto split() const {
            parts_type parts;
            const auto& product = interval_product();
            for (std::size_t i = 0; i < product.size(); i++)
                parts[i] = std::make_shared<HyperCube>(std::move(product[i]));

            return parts;
        }

        auto volume() const {
            long double vol = 1;
            for (const auto& [a, b] : intervals)
                vol *= b - a;
            return vol;
        }

        std::ostream& write(std::ostream& out) const {
            out << "HyperCube<" << N << ", " << typeid(T).name() << ">( ";
            for (const auto& [a, b]: intervals)
                out << "[" << a << ", " << b << "] ";
            out << ")";
            return out;
        }

        friend std::ostream &operator<<(std::ostream &os, const HyperCube& hc) {
            return hc.write(os);
        }


    protected:
        // Split each interval in half, return cartesian product of the halves
        std::vector<intervals_type> interval_product() const {
            std::size_t pos = 0;
            std::vector<intervals_type> rv {{}};
            for (const auto& [a, b] : intervals) {
                std::vector<intervals_type> r;
                T center = (a + b) / 2;

                for (const auto& x : rv)  {
                    // First half
                    r.push_back(x);
                    r.back()[pos] = {a, center};
                    // Second half
                    r.push_back(x);
                    r.back()[pos] = {center, b};
                }

                rv = std::move(r);
                pos++;
            }

            return rv;
        }
    };
}

#endif // HYPERCUBE_H
