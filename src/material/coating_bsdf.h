#pragma once

#include "bsdf.h"
#include "func/frame.h"
#include "func/coating_layer.h"

#include <memory>
#include <utility>

namespace LinRT {
    // Mitsuba coating（SmoothCoating）：光滑电介质涂层包一层任意 nested BSDF。
    // 分量 0 .. N-1 是内层，最后一分量是涂层外表面的 delta 反射。
    // 内层多次反射被丢掉；漫反射基底应使用 plastic，不要用本类去镀。
    class coating_bsdf : public bsdf {
    public:
        coating_bsdf(global_precision int_ior, global_precision ext_ior, std::shared_ptr<bsdf> nested,
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
            , specular_component(this->nested->component_count()) {}

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

            const bool allow_specular = specular_allowed(ctx);
            const bool allow_nested = nested_allowed(ctx);
            if(!allow_specular && !allow_nested) {
                return color(0.0f);
            }

            const auto incoming = coating_refract_in(sample_record.wi, eta, inv_eta);
            const auto prob_specular = coating_specular_probability(incoming.reflectance, specular_sampling_weight);

            auto lobe_sample = sample;
            auto chose_specular = allow_specular;
            if(allow_specular && allow_nested) {
                if(lobe_sample.x < prob_specular) {
                    lobe_sample.x /= prob_specular;
                } else {
                    lobe_sample.x = (lobe_sample.x - prob_specular) / (1.0f - prob_specular);
                    chose_specular = false;
                }
            }

            if(chose_specular) {
                sample_record.wo = coating_reflect(sample_record.wi);
                sample_record.sampled_type = specular_lobe;
                sample_record.sampled_component = specular_component;
                sample_record.pdf = allow_nested ? prob_specular : 1.0f;
                return specular_reflectance * (incoming.reflectance / sample_record.pdf);
            }

            if(incoming.reflectance == 1.0f) {
                return color(0.0f);
            }

            const auto wi_outer = sample_record.wi;
            sample_record.wi = incoming.direction;
            auto weight = nested->sample(ctx, sample_record, lobe_sample);
            const auto wo_inner = sample_record.wo;
            sample_record.wi = wi_outer;
            if(weight.is_zero()) {
                return color(0.0f);
            }

            const auto outgoing = coating_refract_out(wo_inner, eta, inv_eta);
            sample_record.wo = outgoing.direction;
            if(outgoing.reflectance == 1.0f) {
                return color(0.0f);
            }

            weight *= coating_absorption(sigma_a, thickness, incoming.direction, wo_inner);
            weight *= (1.0f - incoming.reflectance) * (1.0f - outgoing.reflectance);

            if(allow_specular) {
                sample_record.pdf *= 1.0f - prob_specular;
                weight /= 1.0f - prob_specular;
            }

            // 权重里的雅可比与 pdf 相消，只改 pdf。离散内层不乘立体角压缩。
            if(measure_of(sample_record.sampled_type) == measure::solid_angle) {
                sample_record.pdf *= coating_solid_angle_factor(inv_eta, sample_record.wo, wo_inner);
            }
            return weight;
        }

        color eval(const bsdf_context& ctx, const bsdf_sample_record& sample_record, measure meas) const override {
            const bool allow_specular = specular_allowed(ctx);
            const bool allow_nested = nested_allowed(ctx);

            if(meas == measure::discrete && allow_specular
                && delta_aligned(coating_reflect(sample_record.wi), sample_record.wo)) {
                return specular_reflectance * fresnel_dielectric_ext(std::abs(sample_record.wi.z), eta);
            }

            if(!allow_nested) {
                return color(0.0f);
            }

            const auto incoming = coating_refract_in(sample_record.wi, eta, inv_eta);
            const auto outgoing = coating_refract_in(sample_record.wo, eta, inv_eta);
            if(incoming.reflectance == 1.0f || outgoing.reflectance == 1.0f) {
                return color(0.0f);
            }

            auto inner = sample_record;
            inner.wi = incoming.direction;
            inner.wo = outgoing.direction;
            auto result = nested->eval(ctx, inner, meas)
                * ((1.0f - incoming.reflectance) * (1.0f - outgoing.reflectance));
            result *= coating_absorption(sigma_a, thickness, incoming.direction, outgoing.direction);
            if(meas == measure::solid_angle) {
                result *= coating_solid_angle_factor(inv_eta, sample_record.wo, outgoing.direction);
            }
            return result;
        }

        global_precision pdf(const bsdf_context& ctx, const bsdf_sample_record& sample_record, measure meas) const override {
            const bool allow_specular = specular_allowed(ctx);
            const bool allow_nested = nested_allowed(ctx);

            const auto incoming = coating_refract_in(sample_record.wi, eta, inv_eta);
            if(meas == measure::discrete && allow_specular
                && delta_aligned(coating_reflect(sample_record.wi), sample_record.wo)) {
                return allow_nested
                    ? coating_specular_probability(incoming.reflectance, specular_sampling_weight)
                    : 1.0f;
            }

            if(!allow_nested) {
                return 0.0f;
            }

            const auto outgoing = coating_refract_in(sample_record.wo, eta, inv_eta);
            if(incoming.reflectance == 1.0f || outgoing.reflectance == 1.0f) {
                return 0.0f;
            }

            auto inner = sample_record;
            inner.wi = incoming.direction;
            inner.wo = outgoing.direction;
            auto density = nested->pdf(ctx, inner, meas);
            if(meas == measure::solid_angle) {
                density *= coating_solid_angle_factor(inv_eta, sample_record.wo, outgoing.direction);
            }
            if(allow_specular) {
                density *= 1.0f - coating_specular_probability(
                    incoming.reflectance, specular_sampling_weight);
            }
            return density;
        }

    private:
        static constexpr unsigned int specular_lobe = SPECULAR | REFLECTIVE | DIELECTRIC;

        bool specular_allowed(const bsdf_context& ctx) const {
            return (ctx.type_mask & specular_lobe) != 0
                && (ctx.component == -1 || ctx.component == specular_component);
        }

        bool nested_allowed(const bsdf_context& ctx) const {
            return (ctx.type_mask & nested->flags()) != 0
                && (ctx.component == -1 || ctx.component < specular_component);
        }

        std::shared_ptr<bsdf> nested;
        global_precision eta;
        global_precision inv_eta;
        global_precision thickness;
        color sigma_a;
        color specular_reflectance;
        global_precision specular_sampling_weight;
        int specular_component;
    };
}
