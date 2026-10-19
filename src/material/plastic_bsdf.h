#pragma once

#include "bsdf.h"
#include "func/frame.h"
#include "func/fresnel.h"

#include <algorithm>
#include <cmath>

namespace LinRT {
    // Mitsuba 的 plastic（SmoothPlastic）：光滑电介质涂层盖在漫反射基底上。
    // 分量 0 是镜面反射，分量 1 是漫反射。只响应正面，背面返回 0。
    // 射线不进入物体，eta 保持 1。内部来回反射收成下面的闭式，不再逐次追踪。
    class plastic_bsdf : public bsdf {
    public:
        plastic_bsdf(global_precision int_ior, global_precision ext_ior,
                     const color& specular_reflectance = color(1.0f),
                     const color& diffuse_reflectance = color(0.5f),
                     bool nonlinear = false)
            : eta(int_ior / ext_ior)
            , inv_eta2(1.0f / (eta * eta))
            , fdr_int(fresnel_diffuse_reflectance(1.0f / eta))
            , specular_reflectance(specular_reflectance)
            , diffuse_reflectance(diffuse_reflectance)
            , specular_sampling_weight(luminance(specular_reflectance)
                / (luminance(diffuse_reflectance) + luminance(specular_reflectance)))
            , nonlinear(nonlinear) {}

        unsigned int flags() const override {
            return specular_lobe | diffuse_lobe;
        }

        int component_count() const override {
            return 2;
        }

        unsigned int component_flags(int index) const override {
            return index == 0 ? specular_lobe : diffuse_lobe;
        }

        color sample(const bsdf_context& ctx, bsdf_sample_record& sample_record, const vec2& sample) const override {
            sample_record.pdf = 0.0f;
            sample_record.eta = 1.0f;
            sample_record.sampled_type = NONE;
            sample_record.sampled_component = -1;

            const bool allow_specular = specular_allowed(ctx);
            const bool allow_diffuse = diffuse_allowed(ctx);
            if((!allow_specular && !allow_diffuse) || sample_record.wi.z <= 0.0f) {
                return color(0.0f);
            }

            const auto fresnel_i = fresnel_dielectric_ext(sample_record.wi.z, eta);
            if(allow_specular && allow_diffuse) {
                // 按菲涅尔反射率抽镜面，再用两个反射率的亮度比把样本推向更亮的那一瓣。
                const auto prob_specular = lobe_probability(fresnel_i);
                if(sample.x < prob_specular) {
                    sample_record.wo = reflect(sample_record.wi, vec3(0.0f, 0.0f, 1.0f));
                    sample_record.sampled_type = specular_lobe;
                    sample_record.sampled_component = 0;
                    sample_record.pdf = prob_specular;
                    return specular_reflectance * (fresnel_i / prob_specular);
                }

                global_precision diffuse_pdf = 0.0f;
                sample_record.wo = cosine_hemisphere(vec2(
                    (sample.x - prob_specular) / (1.0f - prob_specular),
                    sample.y), diffuse_pdf);
                sample_record.pdf = (1.0f - prob_specular) * diffuse_pdf;
                sample_record.sampled_type = diffuse_lobe;
                sample_record.sampled_component = 1;
                const auto fresnel_o = fresnel_dielectric_ext(sample_record.wo.z, eta);
                return diffuse_weight(fresnel_i, fresnel_o) / (1.0f - prob_specular);
            }

            if(allow_specular) {
                sample_record.wo = reflect(sample_record.wi, vec3(0.0f, 0.0f, 1.0f));
                sample_record.sampled_type = specular_lobe;
                sample_record.sampled_component = 0;
                sample_record.pdf = 1.0f;
                return specular_reflectance * fresnel_i;
            }

            global_precision diffuse_pdf = 0.0f;
            sample_record.wo = cosine_hemisphere(sample, diffuse_pdf);
            sample_record.pdf = diffuse_pdf;
            sample_record.sampled_type = diffuse_lobe;
            sample_record.sampled_component = 1;
            const auto fresnel_o = fresnel_dielectric_ext(sample_record.wo.z, eta);
            return diffuse_weight(fresnel_i, fresnel_o);
        }

        color eval(const bsdf_context& ctx, const bsdf_sample_record& sample_record, measure meas) const override {
            const bool allow_specular = specular_allowed(ctx);
            const bool allow_diffuse = diffuse_allowed(ctx);
            if(sample_record.wi.z <= 0.0f) {
                return color(0.0f);
            }

            if(meas == measure::discrete) {
                if(!allow_specular
                    || !delta_aligned(reflect(sample_record.wi, vec3(0.0f, 0.0f, 1.0f)), sample_record.wo)) {
                    return color(0.0f);
                }
                return specular_reflectance * fresnel_dielectric_ext(sample_record.wi.z, eta);
            }

            if(meas != measure::solid_angle || !allow_diffuse || sample_record.wo.z <= 0.0f) {
                return color(0.0f);
            }
            const auto fresnel_i = fresnel_dielectric_ext(sample_record.wi.z, eta);
            const auto fresnel_o = fresnel_dielectric_ext(sample_record.wo.z, eta);
            return diffuse_weight(fresnel_i, fresnel_o) * (sample_record.wo.z / M_PI<global_precision>);
        }

        global_precision pdf(const bsdf_context& ctx, const bsdf_sample_record& sample_record, measure meas) const override {
            const bool allow_specular = specular_allowed(ctx);
            const bool allow_diffuse = diffuse_allowed(ctx);
            if(sample_record.wi.z <= 0.0f) {
                return 0.0f;
            }

            if(meas == measure::discrete) {
                if(!allow_specular
                    || !delta_aligned(reflect(sample_record.wi, vec3(0.0f, 0.0f, 1.0f)), sample_record.wo)) {
                    return 0.0f;
                }
                return allow_diffuse
                    ? lobe_probability(fresnel_dielectric_ext(sample_record.wi.z, eta))
                    : 1.0f;
            }

            if(meas != measure::solid_angle || !allow_diffuse || sample_record.wo.z <= 0.0f) {
                return 0.0f;
            }
            auto prob_specular = 0.0f;
            if(allow_specular) {
                prob_specular = lobe_probability(fresnel_dielectric_ext(sample_record.wi.z, eta));
            }
            return (sample_record.wo.z / M_PI<global_precision>) * (1.0f - prob_specular);
        }

    private:
        static constexpr unsigned int specular_lobe = SPECULAR | REFLECTIVE | DIELECTRIC;
        static constexpr unsigned int diffuse_lobe = DIFFUSE | REFLECTIVE;

        bool specular_allowed(const bsdf_context& ctx) const {
            return (ctx.type_mask & specular_lobe) != 0 && (ctx.component == -1 || ctx.component == 0);
        }

        bool diffuse_allowed(const bsdf_context& ctx) const {
            return (ctx.type_mask & diffuse_lobe) != 0 && (ctx.component == -1 || ctx.component == 1);
        }

        // Mitsuba Spectrum::getLuminance。和 scene::luminance 的系数相同。
        static global_precision luminance(const color& value) {
            return 0.212671f * value.r + 0.715160f * value.g + 0.072169f * value.b;
        }

        global_precision lobe_probability(global_precision fresnel_i) const {
            const auto specular = fresnel_i * specular_sampling_weight;
            const auto diffuse = (1.0f - fresnel_i) * (1.0f - specular_sampling_weight);
            return specular / (specular + diffuse);
        }

        // 基底反照率乘上内部多次反射的和。
        // nonlinear 时每个通道单独用 Kd / (1 - Kd * Fdr_int)，颜色会变暗、饱和度升高。
        // 默认按标量 1 / (1 - Fdr_int) 缩放，贴图颜色保持原样。
        color scaled_diffuse() const {
            if(nonlinear) {
                return color(
                    diffuse_reflectance.r / (1.0f - diffuse_reflectance.r * fdr_int),
                    diffuse_reflectance.g / (1.0f - diffuse_reflectance.g * fdr_int),
                    diffuse_reflectance.b / (1.0f - diffuse_reflectance.b * fdr_int));
            }
            return diffuse_reflectance / (1.0f - fdr_int);
        }

        // 进入涂层的透过率、出去的透过率，以及界面两侧立体角压缩 1/eta^2。
        // 余弦和漫反射 pdf 在 sample 里相消，所以这里还没乘 cos/pi。
        color diffuse_weight(global_precision fresnel_i, global_precision fresnel_o) const {
            return scaled_diffuse() * (inv_eta2 * (1.0f - fresnel_i) * (1.0f - fresnel_o));
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
        global_precision fdr_int;
        color specular_reflectance;
        color diffuse_reflectance;
        global_precision specular_sampling_weight;
        bool nonlinear;
    };
}
