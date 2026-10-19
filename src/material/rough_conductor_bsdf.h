#pragma once

#include "bsdf.h"
#include "func/frame.h"
#include "func/microfacet.h"
#include "func/fresnel.h"

#include <algorithm>
#include <utility>

namespace LinRT {
    // alpha <= 0：理想镜面，对应 Mitsuba 的 conductor。否则是 roughconductor。
    class rough_conductor_bsdf : public bsdf {
    public:
        rough_conductor_bsdf(const conductor_ior& ior, global_precision alpha, const color& specular_reflectance = color(1.0f))
            : eta(ior.eta)
            , k(ior.k)
            , alpha(alpha)
            , specular_reflectance(specular_reflectance)
            , delta(alpha <= 0.0f)
            , lobe(delta ? (SPECULAR | REFLECTIVE | METALLIC) : (GLOSSY | REFLECTIVE | METALLIC))
            , distribution(std::max(alpha, global_precision{0.0001f})) {}

        unsigned int flags() const override {
            return lobe;
        }

        color sample(const bsdf_context& ctx, bsdf_sample_record& sample_record, const vec2& sample) const override {
            sample_record.pdf = 0.0f;
            sample_record.eta = 1.0f; // 反射不改变折射率，俄罗斯轮盘里的 eta 保持原值
            sample_record.sampled_type = NONE;
            sample_record.sampled_component = -1;

            if((ctx.type_mask & lobe) == 0) { return color(0.0f); }
            if(ctx.component != -1 && ctx.component != 0) { return color(0.0f); }
            if(sample_record.wi.z < 0.0f) { return color(0.0f); }

            if(delta) {
                // 几何法线是局部 +Z，反射只翻切向分量。
                sample_record.wo = vec3(-sample_record.wi.x, -sample_record.wi.y, sample_record.wi.z);
                sample_record.sampled_type = lobe;
                sample_record.sampled_component = 0;
                sample_record.pdf = 1.0f;
                return fresnel_conductor(sample_record.wi.z, eta, k) * specular_reflectance;
            }

            const auto m = distribution.sample(sample_record.wi, sample);
            const auto wi_dot_m = dot(sample_record.wi, m);
            if(m.z <= 0.0f || wi_dot_m <= 0.0f) {
                return color(0.0f);
            }

            sample_record.wo = reflect(sample_record.wi, m);
            if(sample_record.wo.z <= 0.0f) {
                return color(0.0f);
            }

            const auto lobe_value = eval_and_pdf(sample_record.wi, sample_record.wo);
            if(lobe_value.second <= 0.0f) {
                return color(0.0f);
            }

            sample_record.pdf = lobe_value.second;
            sample_record.sampled_type = lobe;
            sample_record.sampled_component = 0;
            // 可分离 Smith 加上可见法线抽样时，这个比值等于 F * G1(wo)。
            return lobe_value.first / lobe_value.second;
        }

        color eval(const bsdf_context& ctx, const bsdf_sample_record& sample_record, measure meas) const override {
            if((ctx.type_mask & lobe) == 0) { return color(0.0f); }
            if(ctx.component != -1 && ctx.component != 0) { return color(0.0f); }
            if(delta) {
                if(meas != measure::discrete || sample_record.wo.z <= 0.0f
                    || !delta_aligned(vec3(-sample_record.wi.x, -sample_record.wi.y, sample_record.wi.z),
                                      sample_record.wo)) {
                    return color(0.0f);
                }
                return fresnel_conductor(sample_record.wi.z, eta, k) * specular_reflectance;
            }
            if(meas != measure::solid_angle) { return color(0.0f); }
            return eval_and_pdf(sample_record.wi, sample_record.wo).first;
        }

        global_precision pdf(const bsdf_context& ctx, const bsdf_sample_record& sample_record, measure meas) const override {
            if((ctx.type_mask & lobe) == 0) { return 0.0f; }
            if(ctx.component != -1 && ctx.component != 0) { return 0.0f; }
            if(delta) {
                if(meas != measure::discrete || sample_record.wo.z <= 0.0f
                    || !delta_aligned(vec3(-sample_record.wi.x, -sample_record.wi.y, sample_record.wi.z),
                                      sample_record.wo)) {
                    return 0.0f;
                }
                return 1.0f;
            }
            if(meas != measure::solid_angle) { return 0.0f; }
            return eval_and_pdf(sample_record.wi, sample_record.wo).second;
        }

    private:
        // first 是 f * cosθ_o，second 是 wo 的立体角密度。
        // f = F * D * G / (4 * cosθ_i * cosθ_o)，乘上 cosθ_o 后分母只剩 4 * cosθ_i。
        // 反射时 wo·m = wi·m，从法线密度变到方向密度要除以 4 * wi·m。
        std::pair<color, global_precision> eval_and_pdf(const vec3& wi, const vec3& wo) const {
            if(wi.z <= 1e-6f || wo.z <= 1e-6f) {
                return std::make_pair(color(0.0f), 0.0f);
            }

            const auto half = wi + wo;
            if(half.length_squared() <= 1e-12f) {
                return {color(0.0f), 0.0f};
            }

            const auto m = normalize(half);
            const auto wi_dot_m = dot(wi, m);
            if(m.z <= 0.0f || wi_dot_m <= 0.0f) {
                return {color(0.0f), 0.0f};
            }

            const auto d = distribution.eval(m);
            const auto g = distribution.G(wi, wo, m);
            const auto fresnel = fresnel_conductor(wi_dot_m, eta, k) * specular_reflectance;
            const auto value = fresnel * (d * g / (4.0f * wi.z));

            const auto pdf_m = distribution.pdf(wi, m);
            const auto pdf_wo = pdf_m / (4.0f * wi_dot_m);
            return {value, pdf_wo};
        }

        color eta;
        color k;
        global_precision alpha;
        color specular_reflectance;
        bool delta;
        unsigned int lobe;
        ggx_distribution distribution;
    };
}