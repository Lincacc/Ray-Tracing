#pragma once

#include "bsdf.h"
#include "func/frame.h"
#include "func/fresnel.h"
#include "func/microfacet.h"

#include <algorithm>
#include <utility>

namespace LinRT {
    class rough_dielectric_bsdf : public bsdf {
    public:
        rough_dielectric_bsdf(global_precision int_ior, global_precision ext_ior, global_precision alpha,
                              const color& specular_reflectance = color(1.0f),
                              const color& specular_transmittance = color(1.0f))
            : int_ior(int_ior)
            , ext_ior(ext_ior)
            , alpha(alpha)
            , specular_reflectance(specular_reflectance)
            , specular_transmittance(specular_transmittance)
            , distribution(alpha) {}

        unsigned int flags() const override {
            return GLOSSY | REFLECTIVE | REFRACTIVE | DIELECTRIC;
        }

        color sample(const bsdf_context& ctx, bsdf_sample_record& sample_record, const vec2& sample) const override {
            sample_record.pdf = 0.0f;
            sample_record.eta = 1.0f;
            sample_record.sampled_type = NONE;
            sample_record.sampled_component = -1;

            if(sample_record.wi.z <= 0.0f) {
                return color(0.0f);
            }

            const bool allow_reflection = reflection_allowed(ctx);
            const bool allow_transmission = transmission_allowed(ctx);
            if(!allow_reflection && !allow_transmission) {
                return color(0.0f);
            }

            const auto eta_i = sample_record.front_face ? ext_ior : int_ior;
            const auto eta_t = sample_record.front_face ? int_ior : ext_ior;
            const auto m = distribution.sample(sample_record.wi, sample);
            const auto wi_dot_m = dot(sample_record.wi, m);
            if(m.z <= 0.0f || wi_dot_m <= 0.0f) {
                return color(0.0f);
            }

            // 菲涅尔用的是微表面法线，不是几何法线。全内反射时 reflectance 为 1。
            const auto fresnel = fresnel_dielectric(wi_dot_m, eta_i, eta_t);
            const bool both = allow_reflection && allow_transmission;
            const bool reflect_path = allow_reflection && (!allow_transmission || ctx.lobe_sample < fresnel.reflectance);
            if(reflect_path) {
                sample_record.wo = reflect(sample_record.wi, m);
                if(sample_record.wo.z <= 0.0f) {
                    return color(0.0f);
                }
                const auto lobe = eval_reflection(ctx, sample_record.wi, sample_record.wo, eta_i, eta_t, both);
                if(lobe.second <= 0.0f) {
                    return color(0.0f);
                }
                sample_record.pdf = lobe.second;
                sample_record.sampled_type = GLOSSY | REFLECTIVE | DIELECTRIC;
                sample_record.sampled_component = 0;
                return specular_reflectance * (lobe.first / lobe.second);
            }

            if(fresnel.cos_theta_t <= 0.0f) {
                return color(0.0f);
            }

            sample_record.wo = refract(sample_record.wi, m, fresnel.cos_theta_t, eta_i, eta_t);
            if(sample_record.wo.z >= 0.0f) {
                return color(0.0f);
            }
            sample_record.eta = eta_i / eta_t;
            const auto lobe = eval_transmission(ctx, sample_record.wi, sample_record.wo, eta_i, eta_t, both);
            if(lobe.second <= 0.0f) {
                return color(0.0f);
            }
            sample_record.pdf = lobe.second;
            sample_record.sampled_type = GLOSSY | REFRACTIVE | DIELECTRIC;
            sample_record.sampled_component = 1;
            return specular_transmittance * (lobe.first / lobe.second);
        }

        color eval(const bsdf_context& ctx, const bsdf_sample_record& sample_record) const override {
            return eval_and_pdf(ctx, sample_record).first;
        }

        global_precision pdf(const bsdf_context& ctx, const bsdf_sample_record& sample_record) const override {
            return eval_and_pdf(ctx, sample_record).second;
        }

    private:
        bool reflection_allowed(const bsdf_context& ctx) const {
            return (ctx.type_mask & (GLOSSY | REFLECTIVE | DIELECTRIC)) != 0
                && (ctx.component == -1 || ctx.component == 0);
        }

        bool transmission_allowed(const bsdf_context& ctx) const {
            return (ctx.type_mask & (GLOSSY | REFRACTIVE | DIELECTRIC)) != 0
                && (ctx.component == -1 || ctx.component == 1);
        }

        // first 是 f*cos，second 是 sample() 抽到这个 wo 的立体角密度，含波瓣概率。
        std::pair<color, global_precision> eval_and_pdf(const bsdf_context& ctx, const bsdf_sample_record& sample_record) const {
            if(sample_record.wi.z <= 1e-6f || std::abs(sample_record.wo.z) <= 1e-6f) {
                return {color(0.0f), 0.0f};
            }
            const auto eta_i = sample_record.front_face ? ext_ior : int_ior;
            const auto eta_t = sample_record.front_face ? int_ior : ext_ior;
            const bool both = reflection_allowed(ctx) && transmission_allowed(ctx);
            if(sample_record.wi.z * sample_record.wo.z > 0.0f) {
                return eval_reflection(ctx, sample_record.wi, sample_record.wo, eta_i, eta_t, both);
            }
            return eval_transmission(ctx, sample_record.wi, sample_record.wo, eta_i, eta_t, both);
        }

        std::pair<color, global_precision> eval_reflection(const bsdf_context& ctx, const vec3& wi, const vec3& wo,
                                                           global_precision eta_i, global_precision eta_t, bool both) const {
            if(!reflection_allowed(ctx) || wi.z <= 1e-6f || wo.z <= 1e-6f) {
                return {color(0.0f), 0.0f};
            }
            auto m = wi + wo;
            if(m.length_squared() <= 1e-12f) {
                return {color(0.0f), 0.0f};
            }
            m = normalize(m);
            if(m.z < 0.0f) {
                m = -m;
            }
            const auto wi_dot_m = dot(wi, m);
            if(wi_dot_m <= 0.0f) {
                return {color(0.0f), 0.0f};
            }

            const auto fresnel = fresnel_dielectric(wi_dot_m, eta_i, eta_t).reflectance;
            const auto d = distribution.eval(m);
            const auto g = distribution.G(wi, wo, m);
            const auto pdf_m = distribution.pdf(wi, m);
            // 两个波瓣都开时，抽中反射的概率是 F，密度里要乘 F。只开反射时概率是 1。
            const auto lobe = both ? fresnel : 1.0f;
            const auto value = fresnel * d * g / (4.0f * wi.z);
            const auto pdf_wo = pdf_m * lobe / (4.0f * wi_dot_m);
            return {value, pdf_wo};
        }

        std::pair<color, global_precision> eval_transmission(const bsdf_context& ctx, const vec3& wi, const vec3& wo,
                                                             global_precision eta_i, global_precision eta_t, bool both) const {
            if(!transmission_allowed(ctx) || wi.z <= 1e-6f || wo.z >= -1e-6f || eta_t <= 0.0f) {
                return {color(0.0f), 0.0f};
            }

            // refract 推出 m 与 -(η_i wi + η_t wo) 平行。翻到 z > 0，和可见法线抽样同一侧。
            auto m = -(eta_i * wi + eta_t * wo);
            if(m.length_squared() <= 1e-12f) {
                return {color(0.0f), 0.0f};
            }
            m = normalize(m);
            if(m.z < 0.0f) {
                m = -m;
            }
            const auto wi_dot_m = dot(wi, m);
            const auto wo_dot_m = dot(wo, m);
            if(wi_dot_m <= 0.0f) {
                return {color(0.0f), 0.0f};
            }

            const auto fresnel = fresnel_dielectric(wi_dot_m, eta_i, eta_t).reflectance;
            const auto denom = eta_i * wi_dot_m + eta_t * wo_dot_m;
            if(fresnel >= 1.0f || denom * denom <= 1e-12f) {
                return {color(0.0f), 0.0f};
            }

            const auto d = distribution.eval(m);
            const auto g = distribution.G(wi, wo, m);
            const auto pdf_m = distribution.pdf(wi, m);
            const auto abs_wo_m = std::abs(wo_dot_m);
            // RADIANCE 用 η_i²，和光滑玻璃的 (η_i/η_t)² 一致。IMPORTANCE 用 η_t²，该项在权重里约掉。
            const auto eta_scale = ctx.mode == bsdf_context::transport_mode::RADIANCE ? eta_i : eta_t;
            const auto lobe = both ? (1.0f - fresnel) : 1.0f;
            const auto value = (1.0f - fresnel) * d * g * wi_dot_m * abs_wo_m * eta_scale * eta_scale
                / (wi.z * denom * denom);
            // 半程向量到折射方向的雅可比：η_t² |wo·m| / (η_i wi·m + η_t wo·m)²
            const auto pdf_wo = pdf_m * lobe * (eta_t * eta_t) * abs_wo_m / (denom * denom);
            return {value, pdf_wo};
        }

        global_precision int_ior;
        global_precision ext_ior;
        global_precision alpha;
        color specular_reflectance;
        color specular_transmittance;
        ggx_distribution distribution;
    };
}