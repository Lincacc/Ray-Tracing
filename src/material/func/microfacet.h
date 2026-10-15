#pragma once

#include "core/color.h"

#include <algorithm>
#include <cmath>

namespace LinRT {
    inline vec3 reflect(const vec3& w, const vec3& n) {
        return 2.0f * dot(w, n) * n - w;
    }

    class ggx_distribution {
    public:
        explicit ggx_distribution(global_precision alpha)
        : alpha(std::max(alpha, global_precision{0.0001f})) {}

        // 法线分布 D(m)。朝下的微表面不参与反射。
        global_precision eval(const vec3& m) const {
            if(m.z <= 0.0f) {
                return 0.0f;
            }

            const auto alpha2 = alpha * alpha;
            const auto temp = (m.x * m.x + m.y * m.y) / alpha2 + m.z * m.z;
            return 1.0f / (M_PI<global_precision> * alpha2 * temp * temp);
        }

        // Smith 单方向遮蔽。v 与 m 朝向相反时，这条微表面从 v 看不见。
        global_precision smith_g1(const vec3& v, const vec3& m) const {
            if(v.z < 1e-6f || dot(v, m) * v.z <= 0.0f) {
                return 0.0f;
            }

            const auto tan_theta_alpha_2 = (alpha * alpha) * (v.x * v.x + v.y * v.y) / (v.z * v.z);
            const auto lambda = 0.5f * (-1.0f + std::sqrt(1.0f + tan_theta_alpha_2));
            return 1.0f / (1.0f + lambda);
        }

        // Mitsuba 用可分离 Smith：G = G1(wi) * G1(wo)。
        global_precision G(const vec3& wi, const vec3& wo, const vec3& m) const {
            return smith_g1(wi, m) * smith_g1(wo, m);
        }

        // 可见法线密度：G1(wi) * D(m) * max(wi·m, 0) / cosθ_i
        global_precision pdf(const vec3& wi, const vec3& m) const {
            if(wi.z <= 1e-6f) {
                return 0.0f;
            }
            
            return smith_g1(wi, m) * std::max(0.0f, dot(wi, m)) * eval(m) / wi.z;
        }

        // Heitz 2018：只抽 wi 看得见的微表面法线。wi 和返回值都在局部空间。
        vec3 sample(const vec3& wi, const vec2& u) const {
            // alpha 越小，观察方向被拉得越扁，抽到的法线越靠近镜面反射方向。
            const auto vh = normalize(vec3(alpha * wi.x, alpha * wi.y, wi.z));

            const auto lensq = vh.x * vh.x + vh.y * vh.y;
            const auto t1_axis = lensq > 0.0f
                ? vec3(-vh.y, vh.x, 0.0f) / std::sqrt(lensq)
                : vec3(1.0f, 0.0f, 0.0f);
                const auto t2_axis = cross(vh, t1_axis);

            const auto r = std::sqrt(u.x);
            const auto phi = 2.0f * M_PI<global_precision> * u.y;
            const auto disk_x = r * std::cos(phi);
            auto disk_y = r * std::sin(phi);
            const auto stretch = 0.5f * (1.0f + vh.z);
            disk_y = (1.0f - stretch) * std::sqrt(std::max(0.0f, 1.0f - disk_x * disk_x)) + stretch * disk_y;

            const auto nh = disk_x * t1_axis
                + disk_y * t2_axis
                + std::sqrt(std::max(0.0f, 1.0f - disk_x * disk_x - disk_y * disk_y)) * vh;

            // nh.z 因舍入可能为负，夹到 0，避免法线翻到背面。
            return normalize(vec3(alpha * nh.x, alpha * nh.y, std::max(0.0f, nh.z)));
        }
    private:
        global_precision alpha;
    };
}