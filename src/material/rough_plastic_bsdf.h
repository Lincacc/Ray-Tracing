#pragma once

#include "bsdf.h"
#include "func/frame.h"
#include "func/fresnel.h"
#include "func/microfacet.h"
#include "func/rough_transmittance.h"

#include <algorithm>
#include <cmath>
#include <utility>

namespace LinRT {
    // Mitsuba roughplastic。分量 0 是 GGX 光泽反射，分量 1 是粗糙涂层下的漫反射。
    // 只响应正面。射线留在外面，eta 保持 1。
    class rough_plastic_bsdf : public bsdf {
    public:
        rough_plastic_bsdf(global_precision int_ior, global_precision ext_ior, global_precision alpha,
                           const color& specular_reflectance = color(1.0f),
                           const color& diffuse_reflectance = color(0.5f),
                           bool nonlinear = false)
            : eta(int_ior / ext_ior)
            , inv_eta2(1.0f / (eta * eta))
            , alpha(alpha)
            , specular_reflectance(specular_reflectance)
            , diffuse_reflectance(diffuse_reflectance)
            , specular_sampling_weight(luminance(specular_reflectance)
                / (luminance(diffuse_reflectance) + luminance(specular_reflectance)))
            , nonlinear(nonlinear)
            , distribution(alpha)
            , external(eta, alpha)       // RoughTransmittance::setEta(m_eta)，alpha 是常数
            , internal(1.0f / eta, alpha) // setEta(1/m_eta)，只给内部漫反射 Fdr 用
        {}

        unsigned int flags() const override {
            return glossy_lobe | diffuse_lobe;
        }

        int component_count() const override { return 2; }

        unsigned int component_flags(int index) const override {
            return index == 0 ? glossy_lobe : diffuse_lobe;
        }

        color sample(const bsdf_context& ctx, bsdf_sample_record& sample_record, const vec2& sample) const override {
            sample_record.pdf = 0.0f;
            sample_record.eta = 1.0f;
            sample_record.sampled_type = NONE;
            sample_record.sampled_component = -1;

            const bool allow_glossy = glossy_allowed(ctx);
            const bool allow_diffuse = diffuse_allowed(ctx);
            if((!allow_glossy && !allow_diffuse) || sample_record.wi.z <= 0.0f) {
                return color(0.0f);
            }

            // 光滑 plastic 用 sample.x 选瓣。roughplastic.cpp 用 sample.y，选完再把 y 映射回 [0,1)。
            auto lobe_sample = sample;
            auto chose_glossy = allow_glossy;
            if(allow_glossy && allow_diffuse) {
                const auto prob_glossy = glossy_probability(sample_record.wi.z);
                if(lobe_sample.y < prob_glossy) {
                    lobe_sample.y /= prob_glossy;
                } else {
                    lobe_sample.y = (lobe_sample.y - prob_glossy) / (1.0f - prob_glossy);
                    chose_glossy = false;
                }
            }

            if(chose_glossy) {
                const auto m = distribution.sample(sample_record.wi, lobe_sample);
                sample_record.wo = reflect(sample_record.wi, m);
                sample_record.sampled_type = glossy_lobe;
                sample_record.sampled_component = 0;
                if(sample_record.wo.z <= 0.0f) {
                    return color(0.0f);
                }
            } else {
                global_precision cosine_pdf = 0.0f;
                sample_record.wo = cosine_hemisphere(lobe_sample, cosine_pdf);
                sample_record.sampled_type = diffuse_lobe;
                sample_record.sampled_component = 1;
            }

            // 两瓣都开着时，eval 和 pdf 都含两瓣。权重是二者之和除以混合密度。
            const auto density = pdf(ctx, sample_record);
            if(density == 0.0f) {
                return color(0.0f);
            }
            sample_record.pdf = density;
            return eval(ctx, sample_record) / density;
        }

        color eval(const bsdf_context& ctx, const bsdf_sample_record& sample_record) const override {
            const bool allow_glossy = glossy_allowed(ctx);
            const bool allow_diffuse = diffuse_allowed(ctx);
            if(sample_record.wi.z <= 0.0f || sample_record.wo.z <= 0.0f || (!allow_glossy && !allow_diffuse)) {
                return color(0.0f);
            }

            color result(0.0f);
            if(allow_glossy) {
                const auto half = normalize(sample_record.wi + sample_record.wo);
                const auto d = distribution.eval(half);
                const auto fresnel = fresnel_dielectric_ext(dot(sample_record.wi, half), eta);
                const auto g = distribution.G(sample_record.wi, sample_record.wo, half);
                // f·cos = F D G / (4 cosθi)。分母里的 cosθo 已经和反射 BRDF 约掉。
                result += specular_reflectance * (fresnel * d * g / (4.0f * sample_record.wi.z));
            }

            if(allow_diffuse) {
                // T12、T21 都查外部表（eta = int/ext）。Fdr 才用内部表。
                const auto t12 = external.eval(sample_record.wi.z);
                const auto t21 = external.eval(sample_record.wo.z);
                const auto fdr = 1.0f - internal.eval_diffuse();
                result += scaled_diffuse(fdr) * ((sample_record.wo.z / M_PI<global_precision>) * t12 * t21 * inv_eta2);
            }
            return result;
        }

        global_precision pdf(const bsdf_context& ctx, const bsdf_sample_record& sample_record) const override {
            const bool allow_glossy = glossy_allowed(ctx);
            const bool allow_diffuse = diffuse_allowed(ctx);
            if(sample_record.wi.z <= 0.0f || sample_record.wo.z <= 0.0f || (!allow_glossy && !allow_diffuse)) {
                return 0.0f;
            }

            auto prob_glossy = 1.0f;
            auto prob_diffuse = 1.0f;
            if(allow_glossy && allow_diffuse) {
                prob_glossy = glossy_probability(sample_record.wi.z);
                prob_diffuse = 1.0f - prob_glossy;
            }

            auto density = 0.0f;
            if(allow_glossy) {
                const auto half = normalize(sample_record.wi + sample_record.wo);
                const auto dwh_dwo = 1.0f / (4.0f * dot(sample_record.wo, half));
                density += distribution.pdf(sample_record.wi, half) * dwh_dwo * prob_glossy;
            }
            if(allow_diffuse) {
                density += prob_diffuse * (sample_record.wo.z / M_PI<global_precision>);
            }
            return density;
        }

    private:
        static constexpr unsigned int glossy_lobe = GLOSSY | REFLECTIVE | DIELECTRIC;
        static constexpr unsigned int diffuse_lobe = DIFFUSE | REFLECTIVE;

        bool glossy_allowed(const bsdf_context& ctx) const {
            return (ctx.type_mask & glossy_lobe) != 0 && (ctx.component == -1 || ctx.component == 0);
        }

        bool diffuse_allowed(const bsdf_context& ctx) const {
            return (ctx.type_mask & diffuse_lobe) != 0 && (ctx.component == -1 || ctx.component == 1);
        }

        static global_precision luminance(const color& value) {
            return 0.212671f * value.r + 0.715160f * value.g + 0.072169f * value.b;
        }

        // 1-T 是该入射角的方向半球反射率，再按两个反射率的亮度比重分配样本。
        global_precision glossy_probability(global_precision cos_theta_i) const {
            const auto reflectance = 1.0f - external.eval(cos_theta_i);
            const auto glossy = reflectance * specular_sampling_weight;
            const auto diffuse = (1.0f - reflectance) * (1.0f - specular_sampling_weight);
            return glossy / (glossy + diffuse);
        }

        color scaled_diffuse(global_precision fdr) const {
            if(nonlinear) {
                return color(
                    diffuse_reflectance.r / (1.0f - diffuse_reflectance.r * fdr),
                    diffuse_reflectance.g / (1.0f - diffuse_reflectance.g * fdr),
                    diffuse_reflectance.b / (1.0f - diffuse_reflectance.b * fdr));
            }
            return diffuse_reflectance / (1.0f - fdr);
        }

        static vec3 cosine_hemisphere(const vec2& sample, global_precision& pdf) {
            const auto cos_theta = std::sqrt(sample.x);
            const auto sin_theta = std::sqrt(1.0f - sample.x);
            const auto phi = 2.0f * M_PI<global_precision> * sample.y;
            pdf = cos_theta / M_PI<global_precision>;
            return vec3(std::cos(phi) * sin_theta, std::sin(phi) * sin_theta, cos_theta);
        }

        global_precision eta;
        global_precision inv_eta2;
        global_precision alpha;
        color specular_reflectance;
        color diffuse_reflectance;
        global_precision specular_sampling_weight;
        bool nonlinear;
        ggx_distribution distribution;
        rough_transmittance external;
        rough_transmittance internal;
    };
}