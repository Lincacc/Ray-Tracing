#pragma once

#include "core/core.h"

#include <algorithm>
#include <cstddef>

namespace LinRT {
    // Mitsuba src/libcore/spline.cpp。knots 均匀，导数用中心差分。
    inline global_precision eval_cubic_interp_1d(global_precision x, const global_precision* values,
                                                 std::size_t size, global_precision min, global_precision max) {
        if(!(x >= min && x <= max)) {
            return 0.0f;
        }

        auto t = ((x - min) * static_cast<global_precision>(size - 1)) / (max - min);
        auto k = static_cast<std::size_t>(t);
        if(k > size - 2) {
            k = size - 2;
        }

        const auto f0 = values[k];
        const auto f1 = values[k + 1];
        const auto d0 = k > 0 ? 0.5f * (values[k + 1] - values[k - 1]) : (values[k + 1] - values[k]);
        const auto d1 = (k + 2 < size) ? 0.5f * (values[k + 2] - values[k]) : (values[k + 1] - values[k]);

        t -= static_cast<global_precision>(k);
        const auto t2 = t * t;
        const auto t3 = t2 * t;
        return (2.0f * t3 - 3.0f * t2 + 1.0f) * f0
            + (-2.0f * t3 + 3.0f * t2) * f1
            + (t3 - 2.0f * t2 + t) * d0
            + (t3 - t2) * d1;
    }

    inline global_precision eval_cubic_interp_2d(const global_precision p[2], const global_precision* values,
                                                 const std::size_t size[2],
                                                 const global_precision min[2], const global_precision max[2]) {
        global_precision knot_weights[2][4];
        std::size_t knot[2];

        for(int dim = 0; dim < 2; ++dim) {
            auto* weights = knot_weights[dim];
            if(!(p[dim] >= min[dim] && p[dim] <= max[dim])) {
                return 0.0f;
            }

            auto t = ((p[dim] - min[dim]) * static_cast<global_precision>(size[dim] - 1))
                / (max[dim] - min[dim]);
            knot[dim] = std::min(static_cast<std::size_t>(t), size[dim] - 2);
            t -= static_cast<global_precision>(knot[dim]);

            const auto t2 = t * t;
            const auto t3 = t2 * t;
            weights[0] = 0.0f;
            weights[1] = 2.0f * t3 - 3.0f * t2 + 1.0f;
            weights[2] = -2.0f * t3 + 3.0f * t2;
            weights[3] = 0.0f;
            const auto d0 = t3 - 2.0f * t2 + t;
            const auto d1 = t3 - t2;

            if(knot[dim] > 0) {
                weights[2] += 0.5f * d0;
                weights[0] -= 0.5f * d0;
            } else {
                weights[2] += d0;
                weights[1] -= d0;
            }
            if(knot[dim] + 2 < size[dim]) {
                weights[3] += 0.5f * d1;
                weights[1] -= 0.5f * d1;
            } else {
                weights[2] += d1;
                weights[1] -= d1;
            }
        }

        auto result = 0.0f;
        for(int y = -1; y <= 2; ++y) {
            const auto wy = knot_weights[1][y + 1];
            for(int x = -1; x <= 2; ++x) {
                const auto wxy = knot_weights[0][x + 1] * wy;
                if(wxy == 0.0f) {
                    continue;
                }
                const auto pos = (knot[1] + static_cast<std::size_t>(y)) * size[0]
                    + knot[0] + static_cast<std::size_t>(x);
                result += values[pos] * wxy;
            }
        }
        return result;
    }

    inline global_precision eval_cubic_interp_3d(const global_precision p[3], const global_precision* values,
                                                 const std::size_t size[3],
                                                 const global_precision min[3], const global_precision max[3]) {
        global_precision knot_weights[3][4];
        std::size_t knot[3];

        for(int dim = 0; dim < 3; ++dim) {
            auto* weights = knot_weights[dim];
            if(!(p[dim] >= min[dim] && p[dim] <= max[dim])) {
                return 0.0f;
            }

            auto t = ((p[dim] - min[dim]) * static_cast<global_precision>(size[dim] - 1))
                / (max[dim] - min[dim]);
            knot[dim] = std::min(static_cast<std::size_t>(t), size[dim] - 2);
            t -= static_cast<global_precision>(knot[dim]);

            const auto t2 = t * t;
            const auto t3 = t2 * t;
            weights[0] = 0.0f;
            weights[1] = 2.0f * t3 - 3.0f * t2 + 1.0f;
            weights[2] = -2.0f * t3 + 3.0f * t2;
            weights[3] = 0.0f;
            const auto d0 = t3 - 2.0f * t2 + t;
            const auto d1 = t3 - t2;

            if(knot[dim] > 0) {
                weights[2] += 0.5f * d0;
                weights[0] -= 0.5f * d0;
            } else {
                weights[2] += d0;
                weights[1] -= d0;
            }
            if(knot[dim] + 2 < size[dim]) {
                weights[3] += 0.5f * d1;
                weights[1] -= 0.5f * d1;
            } else {
                weights[2] += d1;
                weights[1] -= d1;
            }
        }

        auto result = 0.0f;
        for(int z = -1; z <= 2; ++z) {
            const auto wz = knot_weights[2][z + 1];
            for(int y = -1; y <= 2; ++y) {
                const auto wyz = knot_weights[1][y + 1] * wz;
                for(int x = -1; x <= 2; ++x) {
                    const auto wxyz = knot_weights[0][x + 1] * wyz;
                    if(wxyz == 0.0f) {
                        continue;
                    }
                    const auto pos = ((knot[2] + static_cast<std::size_t>(z)) * size[1]
                        + (knot[1] + static_cast<std::size_t>(y))) * size[0]
                        + knot[0] + static_cast<std::size_t>(x);
                    result += values[pos] * wxyz;
                }
            }
        }
        return result;
    }
}