#pragma once

#include "bsdf.h"
#include "func/frame.h"
#include "func/coating_layer.h"
#include "func/microfacet.h"
#include "func/rough_transmittance.h"

#include <memory>
#include <utility>

namespace LinRT {
    // Mitsuba roughcoating（RoughCoating）：粗糙电介质涂层包一层任意 nested BSDF。
    // 分量 0 .. N-1 是内层，最后一分量是涂层外表面的 GGX 光泽反射。
    // 粗糙只作用在外表面反射和界面透过率上；内层仍按光滑涂层做方向变换。
    // 涂层内的多次内反射被丢掉。
    class rough_coating_bsdf : public bsdf {
    public:
        rough_coating_bsdf(global_precision int_ior, global_precision ext_ior, global_precision alpha,
                           std::shared_ptr<bsdf> nested,
                           const color& specular_reflectance = color(1.0f),
                           const color& sigma_a = color(0.0f),
                           global_precision thickness = 1.0f)
            : nested(std::move(nested))
            , eta(int_ior / ext_ior)
            , inv_eta(ext_ior / int_ior)
            , thickness(thickness)
            , sigma_a(sigma_a)
            , specular_reflectance(specular_reflectance)
            , specular_sampling_weight(coating_specular_sampling_weight(sigma_a, thickness))
            , specular_component(this->nested->component_count())
            , distribution(alpha)
            , transmittance(int_ior / ext_ior, alpha) {}

        unsigned int flags() const override {
            return nested->flags() | specular_lobe;
        }

        int component_count() const override {
            return specular_component + 1;
        }

        unsigned int component_flags(int index) const override {
            return index == specular_component ? specular_lobe : nested->component_flags(index);
        }

        color diffuse_reflectance() const override {
            return nested->diffuse_reflectance();
        }

        color sample(const bsdf_context& ctx, bsdf_sample_record& sample_record, const vec2& sample) const override {
            sample_record.pdf = 0.0f;
            sample_record.eta = 1.0f;
            sample_record.sampled_type = NONE;
            sample_record.sampled_component = -1;

            const bool allow_specular = specular_component_allowed(ctx);
            const bool allow_nested = nested_allowed(ctx);
            if(!allow_specular && !allow_nested) {
                return color(0.0f);
            }

            // roughcoating.cpp 用 sample.y 选瓣，与 roughplastic 相同。
            auto lobe_sample = sample;
            auto chose_specular = allow_specular;
            if(allow_specular && allow_nested) {
                const auto prob_specular = specular_probability(sample_record.wi.z);
                if(lobe_sample.y < prob_specular) {
                    lobe_sample.y /= prob_specular;
                } else {
                    lobe_sample.y = (lobe_sample.y - prob_specular) / (1.0f - prob_specular);
                    chose_specular = false;
                }
            }

            if(chose_specular) {
                const auto m = distribution.sample(sample_record.wi, lobe_sample);
                sample_record.wo = reflect(sample_record.wi, m);
                sample_record.sampled_type = specular_lobe;
                sample_record.sampled_component = specular_component;
                sample_record.eta = 1.0f;
                if(sample_record.wo.z * sample_record.wi.z <= 0.0f) {
                    return color(0.0f);
                }
            } else {
                const auto wi_outer = sample_record.wi;
                sample_record.wi = coating_refract_to(
                    coating_destination::interior, wi_outer, eta, inv_eta);
                auto weight = nested->sample(ctx, sample_record, lobe_sample);
                sample_record.wi = wi_outer;
                if(weight.is_zero()) {
                    return color(0.0f);
                }
                sample_record.wo = coating_refract_to(
                    coating_destination::exterior, sample_record.wo, eta, inv_eta);
                if(sample_record.wo.is_zero()) {
                    return color(0.0f);
                }
            }

            // 选瓣之后按 eval/pdf 混两路，与 roughcoating.cpp 的 sample(bRec, pdf, sample) 相同。
            const auto meas = measure_of(sample_record.sampled_type);
            const auto density = pdf(ctx, sample_record, meas);
            if(density == 0.0f) {
                return color(0.0f);
            }
            sample_record.pdf = density;
            return eval(ctx, sample_record, meas) / density;
        }

        color eval(const bsdf_context& ctx, const bsdf_sample_record& sample_record, measure meas) const override {
            const bool allow_nested = nested_allowed(ctx);
            const bool allow_specular = specular_lobe_allowed(ctx, meas);

            color result(0.0f);
            if(allow_specular && sample_record.wo.z * sample_record.wi.z > 0.0f) {
                auto half = sample_record.wi + sample_record.wo;
                if(half.length_squared() > 1e-12f) {
                    half = normalize(half) * signum(sample_record.wo.z);
                    const auto d = distribution.eval(half);
                    const auto fresnel = fresnel_dielectric_ext(std::abs(dot(sample_record.wi, half)), eta);
                    const auto g = distribution.G(sample_record.wi, sample_record.wo, half);
                    // f·cos = F D G / (4 |cosθi|)。背面用绝对值，与 roughcoating.cpp 一致。
                    result += specular_reflectance
                        * (fresnel * d * g / (4.0f * std::abs(sample_record.wi.z)));
                }
            }

            if(allow_nested) {
                result += nested_eval(ctx, sample_record, meas);
            }
            return result;
        }

        global_precision pdf(const bsdf_context& ctx, const bsdf_sample_record& sample_record, measure meas) const override {
            const bool allow_nested = nested_allowed(ctx);
            const bool allow_specular = specular_lobe_allowed(ctx, meas);

            auto prob_specular = 1.0f;
            auto prob_nested = 1.0f;
            if(allow_specular && allow_nested) {
                prob_specular = specular_probability(sample_record.wi.z);
                prob_nested = 1.0f - prob_specular;
            }

            auto density = 0.0f;
            if(allow_specular && sample_record.wo.z * sample_record.wi.z > 0.0f) {
                auto half = sample_record.wi + sample_record.wo;
                if(half.length_squared() > 1e-12f) {
                    half = normalize(half) * signum(sample_record.wo.z);
                    const auto dwh_dwo = 1.0f / (4.0f * std::abs(dot(sample_record.wo, half)));
                    density += distribution.pdf(sample_record.wi, half) * dwh_dwo * prob_specular;
                }
            }

            if(allow_nested) {
                density += nested_pdf(ctx, sample_record, meas) * prob_nested;
            }
            return density;
        }

    private:
        static constexpr unsigned int specular_lobe = GLOSSY | REFLECTIVE | DIELECTRIC;

        bool specular_component_allowed(const bsdf_context& ctx) const {
            return (ctx.type_mask & specular_lobe) != 0
                && (ctx.component == -1 || ctx.component == specular_component);
        }

        bool specular_lobe_allowed(const bsdf_context& ctx, measure meas) const {
            return specular_component_allowed(ctx) && meas == measure::solid_angle;
        }

        bool nested_allowed(const bsdf_context& ctx) const {
            return (ctx.type_mask & nested->flags()) != 0
                && (ctx.component == -1 || ctx.component < specular_component);
        }

        // 1-T 是该入射角的粗糙界面反射率，再按吸收算出的权重把样本分给外层或内层。
        global_precision specular_probability(global_precision cos_theta_i) const {
            return coating_specular_probability(
                1.0f - transmittance.eval(std::abs(cos_theta_i)), specular_sampling_weight);
        }

        color nested_eval(const bsdf_context& ctx, const bsdf_sample_record& sample_record, measure meas) const {
            auto inner = sample_record;
            inner.wi = coating_refract_to(coating_destination::interior, sample_record.wi, eta, inv_eta);
            inner.wo = coating_refract_to(coating_destination::interior, sample_record.wo, eta, inv_eta);
            // refractTo 全内反射时方向为零，内层贡献为 0，不能再除 cosθ_inner。
            if(inner.wi.is_zero() || inner.wo.is_zero()) {
                return color(0.0f);
            }

            auto result = nested->eval(ctx, inner, meas)
                * transmittance.eval(std::abs(sample_record.wi.z))
                * transmittance.eval(std::abs(sample_record.wo.z));
            result *= coating_absorption(sigma_a, thickness, inner.wi, inner.wo);
            if(meas == measure::solid_angle) {
                result *= coating_solid_angle_factor(inv_eta, sample_record.wo, inner.wo);
            }
            return result;
        }

        global_precision nested_pdf(const bsdf_context& ctx, const bsdf_sample_record& sample_record, measure meas) const {
            auto inner = sample_record;
            inner.wi = coating_refract_to(coating_destination::interior, sample_record.wi, eta, inv_eta);
            inner.wo = coating_refract_to(coating_destination::interior, sample_record.wo, eta, inv_eta);
            if(inner.wi.is_zero() || inner.wo.is_zero()) {
                return 0.0f;
            }

            auto density = nested->pdf(ctx, inner, meas);
            if(meas == measure::solid_angle) {
                density *= coating_solid_angle_factor(inv_eta, sample_record.wo, inner.wo);
            }
            return density;
        }

        std::shared_ptr<bsdf> nested;
        global_precision eta;
        global_precision inv_eta;
        global_precision thickness;
        color sigma_a;
        color specular_reflectance;
        global_precision specular_sampling_weight;
        int specular_component;
        ggx_distribution distribution;
        rough_transmittance transmittance;
    };
}
