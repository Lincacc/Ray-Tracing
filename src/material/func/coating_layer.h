#pragma once

#include "core/color.h"
#include "core/vec3.h"
#include "material/func/fresnel.h"

#include <cmath>

namespace LinRT {
    struct coating_refracted {
        vec3 direction;
        global_precision reflectance = 1.0f;
    };

    // SmoothCoating::reflect：局部 +Z 镜面反射。
    inline vec3 coating_reflect(const vec3& wi) {
        return vec3(-wi.x, -wi.y, wi.z);
    }

    inline global_precision signum(global_precision value) {
        return static_cast<global_precision>((value > 0.0f) - (value < 0.0f));
    }

    // SmoothCoating::refractIn。进入涂层后方向仍与 wi 同半球，xy 按 1/eta 压缩。
    // R 是涂层外表面的菲涅尔反射率；全内反射时 R = 1。
    inline coating_refracted coating_refract_in(const vec3& wi, global_precision eta, global_precision inv_eta) {
        coating_refracted result;
        global_precision cos_theta_t = 0.0f;
        result.reflectance = fresnel_dielectric_ext(std::abs(wi.z), cos_theta_t, eta);
        result.direction = vec3(inv_eta * wi.x, inv_eta * wi.y, -signum(wi.z) * cos_theta_t);
        return result;
    }

    // SmoothCoating::refractOut。从涂层回到外部，xy 按 eta 展开。
    inline coating_refracted coating_refract_out(const vec3& wi, global_precision eta, global_precision inv_eta) {
        coating_refracted result;
        global_precision cos_theta_t = 0.0f;
        result.reflectance = fresnel_dielectric_ext(std::abs(wi.z), cos_theta_t, inv_eta);
        result.direction = vec3(eta * wi.x, eta * wi.y, -signum(wi.z) * cos_theta_t);
        return result;
    }

    // RoughCoating::EDestination。interior：外→涂层；exterior：涂层→外。
    enum class coating_destination {
        interior,
        exterior
    };

    // RoughCoating::refractTo。Snell 只改切向尺度，z 符号与 wi 相同，内层仍在正面求值。
    // 全内反射返回零向量，调用方用 is_zero() 判定。
    inline vec3 coating_refract_to(coating_destination dest, const vec3& wi,
                                   global_precision eta, global_precision inv_eta) {
        const auto scale = dest == coating_destination::interior ? inv_eta : eta;
        const auto sin_theta_t2 = scale * scale * (wi.x * wi.x + wi.y * wi.y);
        if(sin_theta_t2 >= 1.0f) {
            return vec3(0.0f);
        }
        const auto cos_theta_t = std::sqrt(1.0f - sin_theta_t2);
        return vec3(scale * wi.x, scale * wi.y, wi.z > 0.0f ? cos_theta_t : -cos_theta_t);
    }

    // 立体角压缩 1/eta^2，以及把内层 eval 里的 cosθ_inner 换成外层 cosθ_outer。
    inline global_precision coating_solid_angle_factor(global_precision inv_eta, const vec3& wo, const vec3& wo_inner) {
        return inv_eta * inv_eta * wo.z / wo_inner.z;
    }

    // configure()：平均吸收 exp(-2 t σ_a) 的 RGB 均值，用来把样本推向镜面或内层。
    inline global_precision coating_specular_sampling_weight(const color& sigma_a, global_precision thickness) {
        const auto optical = sigma_a * (-2.0f * thickness);
        const auto avg_absorption = (std::exp(optical.r) + std::exp(optical.g) + std::exp(optical.b)) / 3.0f;
        return 1.0f / (avg_absorption + 1.0f);
    }

    // 涂层内沿入射、出射各走一段：exp(-σ_a t (1/|cos i| + 1/|cos o|))。
    inline color coating_absorption(const color& sigma_a, global_precision thickness,
                                    const vec3& wi_inner, const vec3& wo_inner) {
        if(sigma_a.is_zero()) {
            return color(1.0f);
        }
        const auto path_length = (1.0f / std::abs(wi_inner.z)) + (1.0f / std::abs(wo_inner.z));
        const auto exponent = sigma_a * (-thickness * path_length);
        return color(std::exp(exponent.r), std::exp(exponent.g), std::exp(exponent.b));
    }

    // 按菲涅尔和外层采样权重重分配样本，与 coating.cpp 的 probSpecular 相同。
    inline global_precision coating_specular_probability(global_precision reflectance,
                                                         global_precision specular_sampling_weight) {
        const auto specular = reflectance * specular_sampling_weight;
        const auto nested = (1.0f - reflectance) * (1.0f - specular_sampling_weight);
        return specular / (specular + nested);
    }
}
