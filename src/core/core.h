#pragma once

#include <limits>

namespace LinRT {
    template <typename T>
    constexpr T M_PI = static_cast<T>(3.14159265358979323846);

    using global_precision = float;

    constexpr global_precision INFINITY = std::numeric_limits<global_precision>::infinity();

    inline global_precision radians(global_precision degrees) {
        return degrees * M_PI<global_precision> / 180.0f;
    }
}