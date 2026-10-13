#pragma once 

#include "core/vec3.h"

namespace LinRT {
    class onb {
    public:
        inline void build_from_normal(const vec3& n) {
            axis[2] = LinRT::normalize(n);

            // 选取一个不平行的辅助向量来构造切线
            vec3 a = (std::abs(axis[2].x) > 0.9f) ? vec3(0.0f, 1.0f, 0.0f) : vec3(1.0f, 0.0f, 0.0f);
            axis[1] = normalize(cross(axis[2], a));
            axis[0] = cross(axis[2], axis[1]);
        }

        inline vec3 local_to_world(const vec3& a) const {
            return a.x * axis[0] + a.y * axis[1] + a.z * axis[2];
        }

        inline vec3 world_to_local(const vec3& a) const {
            return vec3(dot(a, axis[0]), dot(a, axis[1]), dot(a, axis[2]));
        }

    private:
        vec3 axis[3];
    };
}