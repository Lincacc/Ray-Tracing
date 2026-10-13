#pragma once

namespace LinRT {
    template <typename T>
    constexpr T M_PI = static_cast<T>(3.14159265358979323846);

    using global_precision = float;

    inline global_precision radians(global_precision degrees) {
        return degrees * M_PI<global_precision> / 180.0f;
    }
}