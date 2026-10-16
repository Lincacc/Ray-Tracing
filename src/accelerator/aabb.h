#pragma once

#include "core/ray.h"

#include <algorithm>

namespace LinRT {
    class aabb {
    public:
        void extend(const vec3& point) {
            minimum.x = std::min(minimum.x, point.x);
            minimum.y = std::min(minimum.y, point.y);
            minimum.z = std::min(minimum.z, point.z);
            maximum.x = std::max(maximum.x, point.x);
            maximum.y = std::max(maximum.y, point.y);
            maximum.z = std::max(maximum.z, point.z);
        }

        void extend(const aabb& box) {
            extend(box.minimum);
            extend(box.maximum);
        }

        // 与 SurfaceAreaHeuristic3::getQuantity 相同，用表面积。
        global_precision surface_area() const {
            const auto edge = maximum - minimum;
            return 2.0f * (edge.x * edge.y + edge.x * edge.z + edge.y * edge.z);
        }
        
    public:
        vec3 minimum = vec3(INFINITY);
        vec3 maximum = vec3(-INFINITY);
    };
}