#pragma once

#include "bsdf.h"

namespace LinRT {
    class diffuse_bsdf : public bsdf {
    public:
        diffuse_bsdf(const color& albedo) : albedo(albedo) {}

        unsigned int flags() const override {
            return DIFFUSE | REFLECTIVE;
        }

        color sample(const bsdf_context& ctx, bsdf_sample_record& sample_record, const vec2& sample) const override {
            constexpr auto lobe = DIFFUSE | REFLECTIVE;

            sample_record.pdf = 0.0f;
            sample_record.eta = 1.0f;
            sample_record.sampled_type = NONE;
            sample_record.sampled_component = -1;

            if(!(ctx.type_mask & lobe)) { return color(0.0f); }
            if(ctx.component != -1 && ctx.component != 0) { return color(0.0f); }
            if(sample_record.wi.z <= .0f) { return color(0.0f); }

            // 1. 在局部余弦加权半球上采样一个出射方向 wo
            // （利用 sample.x 和 sample.y 进行余弦加权半球采样）
            auto cos_theta = std::sqrt(sample.x);
            auto sin_theta = std::sqrt(1.0f - sample.x);
            auto phi = 2.0f * M_PI<global_precision> * sample.y;

            sample_record.wo = vec3(
                std::cos(phi) * sin_theta, 
                std::sin(phi) * sin_theta, 
                cos_theta);
            sample_record.pdf = cos_theta / M_PI<global_precision>;
            sample_record.sampled_type = lobe;
            sample_record.sampled_component = 0;

            // 2. 返回 f * cos_theta / pdf
            // 对于余弦加权采样：f = albedo / pi, cos_theta = wo.z, pdf = cos_theta / pi
            // 最终约分结果：直接返回 albedo！
            return albedo;
        }

        color eval(const bsdf_context& ctx, const bsdf_sample_record& sample_record) const override {
            constexpr auto lobe = DIFFUSE | REFLECTIVE;

            if(!(ctx.type_mask & lobe)) { return color(0.0f); }
            if(ctx.component != -1 && ctx.component != 0) { return color(0.0f); }
            if(sample_record.wi.z <= .0f || sample_record.wo.z <= .0f) { return color(0.0f); }

            // f_r * cos_theta = (albedo / pi) * cos(theta)
            return albedo * (sample_record.wo.z / M_PI<global_precision>);
        }

        global_precision pdf(const bsdf_context& ctx, const bsdf_sample_record& sample_record) const override {
            constexpr auto lobe = DIFFUSE | REFLECTIVE;

            if(!(ctx.type_mask & lobe)) { return 0.0f; }
            if(ctx.component != -1 && ctx.component != 0) { return 0.0f; }
            if(sample_record.wi.z <= .0f || sample_record.wo.z <= .0f) { return 0.0f; }

            // 余弦加权半球采样的 PDF 为 cos_theta / pi
            return sample_record.wo.z / M_PI<global_precision>;
        }

    private:
        color albedo;
    };
}