#ifndef CONST_H
#define CONST_H

#include <cstdint>
#include <cmath>

namespace Integration {
    constexpr double pi() { return std::atan(1) * 4; }
    enum REGION_STATE : int8_t { REJECTED, CONTAINED, INDEFINITE };
}

#endif // CONST_H
