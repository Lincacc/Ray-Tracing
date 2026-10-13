#pragma once

#include "vec3.h"

namespace LinRT {
    class ray {
    public:
        ray() = default;
        ray(const vec3& origin, const vec3& direction) : origin(origin), direction(direction) {}

        vec3 at(global_precision t) const {
            return origin + t * direction;
        }

    public:
        vec3 origin;
        vec3 direction;
    };
}