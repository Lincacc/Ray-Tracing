#pragma once

#include "core/color.h"
#include "material/func/gauss_lobatto.h"

#include <algorithm>
#include <cmath>
#include <complex>

namespace LinRT {
    struct conductor_ior {
        color eta;
        color k;
        
        static conductor_ior copper() {
            return {
                color(0.200438f, 0.924033f, 1.102212f),
                color(3.912948f, 2.452847f, 2.142188f)
            };
        }

        static conductor_ior gold() {
            return {
                color(0.143119f, 0.374957f, 1.442478f),
                color(3.983160f, 2.385721f, 1.603215f)
            };
        }

        static conductor_ior aluminum() {
            return {
                color(1.345600f, 0.965210f, 0.617220f),
                color(7.474600f, 6.399500f, 5.303100f)
            };
        }
        static conductor_ior silver() {
            return {
                color(0.155265f, 0.116723f, 0.138381f),
                color(4.828350f, 3.122245f, 2.146900f)
            };
        }
    };

    inline global_precision fresnel_conductor(global_precision cos_theta_i, global_precision eta, global_precision k) {
        cos_theta_i = std::clamp(cos_theta_i, 0.0f, 1.0f);

        const std::complex<global_precision> ior(eta, k);
        const auto sin2_theta_i = 1.0f - cos_theta_i * cos_theta_i;
        const auto sin2_theta_t = sin2_theta_i / (ior * ior);
        const auto cos_theta_t = std::sqrt(std::complex<global_precision>(1.0f) - sin2_theta_t);

        const auto r_parallel = (ior * cos_theta_i - cos_theta_t) / (ior * cos_theta_i + cos_theta_t);
        const auto r_perpendicular = (cos_theta_i - ior * cos_theta_t) /(cos_theta_i + ior * cos_theta_t);

        // std::norm 是模的平方，振幅反射率变成能量反射率
        return 0.5f * (std::norm(r_parallel) + std::norm(r_perpendicular));
    }

    inline color fresnel_conductor(global_precision cos_theta_i, const color& eta, const color& k) {
        return color(
            fresnel_conductor(cos_theta_i, eta.r, k.r),
            fresnel_conductor(cos_theta_i, eta.g, k.g),
            fresnel_conductor(cos_theta_i, eta.b, k.b)
        );
    }

    struct dielectric_fresnel_result {
        global_precision reflectance = 1.0f;
        global_precision cos_theta_t = 0.0f;
    };

    // eta_i 是入射侧折射率，eta_t 是折射侧。全内反射时 reflectance 为 1，cos_theta_t 为 0。
    inline dielectric_fresnel_result fresnel_dielectric(global_precision cos_theta_i,
                                                        global_precision eta_i,
                                                        global_precision eta_t) {
        dielectric_fresnel_result result;
        if(eta_i <= 0.0f || eta_t <= 0.0f) {
            return result;
        }

        cos_theta_i = std::clamp(cos_theta_i, 0.0f, 1.0f);
        const auto sin2_i = 1.0f - cos_theta_i * cos_theta_i;
        const auto scale = eta_i / eta_t;
        const auto sin2_t = scale * scale * sin2_i;
        if(sin2_t > 1.0f) {
            return result;
        }

        result.cos_theta_t = std::sqrt(1.0f - sin2_t);
        const auto r_perp = (eta_i * cos_theta_i - eta_t * result.cos_theta_t)
            / (eta_i * cos_theta_i + eta_t * result.cos_theta_t);
        const auto r_parl = (eta_t * cos_theta_i - eta_i * result.cos_theta_t)
            / (eta_t * cos_theta_i + eta_i * result.cos_theta_t);
        result.reflectance = 0.5f * (r_perp * r_perp + r_parl * r_parl);
        return result;
    }

    // Mitsuba fresnelDielectricExt。eta 是相对折射率：cos_theta_i > 0 时，透射侧 / 入射侧。
    // 全内反射返回 1，并把 cos_theta_t 写成 0。eta == 1 时反射率为 0，透射方向反向。
    inline global_precision fresnel_dielectric_ext(global_precision cos_theta_i, global_precision& cos_theta_t,
                                                   global_precision eta) {
        if(eta == 1.0f) {
            cos_theta_t = -cos_theta_i;
            return 0.0f;
        }

        const auto scale = cos_theta_i > 0.0f ? (1.0f / eta) : eta;
        const auto cos_theta_t_sqr = 1.0f - (1.0f - cos_theta_i * cos_theta_i) * (scale * scale);
        if(cos_theta_t_sqr <= 0.0f) {
            cos_theta_t = 0.0f;
            return 1.0f;
        }

        const auto cos_i = std::abs(cos_theta_i);
        const auto cos_t = std::sqrt(cos_theta_t_sqr);
        const auto rs = (cos_i - eta * cos_t) / (cos_i + eta * cos_t);
        const auto rp = (eta * cos_i - cos_t) / (eta * cos_i + cos_t);
        cos_theta_t = cos_theta_i > 0.0f ? -cos_t : cos_t;
        return 0.5f * (rs * rs + rp * rp);
    }

    inline global_precision fresnel_dielectric_ext(global_precision cos_theta_i, global_precision eta) {
        global_precision cos_theta_t = 0.0f;
        return fresnel_dielectric_ext(cos_theta_i, cos_theta_t, eta);
    }

    // 余弦加权半球上的平均菲涅尔反射率。被积函数与积分器参数和 Mitsuba util.cpp 一致：
    // fresnelDielectricExt(sqrt(xi), eta)，xi 从 0 到 1，fast = false。
    inline global_precision fresnel_diffuse_reflectance(global_precision eta) {
        const gauss_lobatto_integrator quadrature(1024, 0.0f, 1e-5f);
        return quadrature.integrate([eta](global_precision xi) {
            return fresnel_dielectric_ext(std::sqrt(xi), eta);
        }, 0.0f, 1.0f);
    }
}