#pragma once

#include "core/vec3.h"

namespace LinRT {
    inline vec3 reflect(const vec3& w, const vec3& n) {
        return 2.0f * dot(w, n) * n - w;
    }

    // wi 和 n 都指向入射侧，cos_theta_t 是透射方向和法线夹角的余弦，取正数。
    // 返回的 wo 指向另一侧。局部法线是 +Z 时，wo.z 为负。
    inline vec3 refract(const vec3& wi, const vec3& n, global_precision cos_theta_t,
                        global_precision eta_i, global_precision eta_t) {
        const auto scale = eta_i / eta_t;
        return n * (scale * dot(wi, n) - cos_theta_t) - wi * scale;
    }
}