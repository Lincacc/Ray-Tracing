#pragma once 

#include "bsdf.h"
#include "func/frame.h"
#include "func/fresnel.h"

namespace LinRT {
    // Mitsuba 的 dielectric。int_ior 是物体内部，ext_ior 是外部，空气用 1
    class dielectric_bsdf : public bsdf {
    public:
        dielectric_bsdf(global_precision int_ior, global_precision ext_ior,
                        const color& specular_reflectance = color(1.0f),
                        const color& specular_transmittance = color(1.0f))
            : int_ior(int_ior)
            , ext_ior(ext_ior)
            , specular_reflectance(specular_reflectance)
            , specular_transmittance(specular_transmittance) {}

        unsigned int flags() const override { 
            return reflection_lobe | transmission_lobe; 
        }

        int component_count() const override {
            return 2;
        }

        unsigned int component_flags(int index) const override {
            return index == 0 ? reflection_lobe : transmission_lobe;
        }

        color sample(const bsdf_context& ctx, bsdf_sample_record& sample_record, const vec2& sample) const override {
            sample_record.pdf = 0.0f;
            sample_record.eta = 1.0f;
            sample_record.sampled_type = NONE;
            sample_record.sampled_component = -1;

            if(sample_record.wi.z <= 0.0f) {
                return color(0.0f);
            }

            const bool allow_reflection = (ctx.type_mask & reflection_lobe) != 0
                && (ctx.component == -1 || ctx.component == 0);
            const bool allow_transmission = (ctx.type_mask & transmission_lobe) != 0
                && (ctx.component == -1 || ctx.component == 1);
            if(!allow_reflection && !allow_transmission) {
                return color(0.0f);
            }

            // 正面：入射侧是外部。背面：射线已经在玻璃里，两侧对调。
            const auto eta_i = sample_record.front_face ? ext_ior : int_ior;
            const auto eta_t = sample_record.front_face ? int_ior : ext_ior;
            const auto fresnel = fresnel_dielectric(sample_record.wi.z, eta_i, eta_t);

            const bool reflect_path = allow_reflection
                && (!allow_transmission || sample.x < fresnel.reflectance);
            if(reflect_path) {
                sample_record.wo = reflect(sample_record.wi, vec3(0.0f, 0.0f, 1.0f));
                sample_record.sampled_type = reflection_lobe;
                sample_record.sampled_component = 0;
                // 离散 pdf：两瓣都开是 F，只开反射是 1。权重 = F·R / pdf。
                sample_record.pdf = allow_transmission ? fresnel.reflectance : 1.0f;
                return allow_transmission ? specular_reflectance
                    : specular_reflectance * fresnel.reflectance;
            }

            if(fresnel.reflectance >= 1.0f) {
                return color(0.0f); // 全内反射，又没允许走反射
            }

            sample_record.wo = refract(sample_record.wi, vec3(0.0f, 0.0f, 1.0f),
                                        fresnel.cos_theta_t, eta_i, eta_t);
            sample_record.eta = eta_i / eta_t;
            sample_record.sampled_type = transmission_lobe;
            sample_record.sampled_component = 1;
            sample_record.pdf = allow_reflection ? (1.0f - fresnel.reflectance) : 1.0f;

            auto weight = specular_transmittance;
            if(!allow_reflection) {
                weight *= 1.0f - fresnel.reflectance;
            }
            if(ctx.mode == bsdf_context::transport_mode::RADIANCE) {
                weight *= sample_record.eta * sample_record.eta;
            }

            return weight;
        }

        color eval(const bsdf_context& ctx, const bsdf_sample_record& sample_record, measure meas) const override {
            if(meas != measure::discrete || sample_record.wi.z <= 0.0f) {
                return color(0.0f);
            }

            const bool allow_reflection = (ctx.type_mask & reflection_lobe) != 0
                && (ctx.component == -1 || ctx.component == 0);
            const bool allow_transmission = (ctx.type_mask & transmission_lobe) != 0
                && (ctx.component == -1 || ctx.component == 1);

            const auto eta_i = sample_record.front_face ? ext_ior : int_ior;
            const auto eta_t = sample_record.front_face ? int_ior : ext_ior;
            const auto fresnel = fresnel_dielectric(sample_record.wi.z, eta_i, eta_t);

            if(sample_record.wi.z * sample_record.wo.z >= 0.0f) {
                if(!allow_reflection
                    || !delta_aligned(reflect(sample_record.wi, vec3(0.0f, 0.0f, 1.0f)), sample_record.wo)) {
                    return color(0.0f);
                }
                return specular_reflectance * fresnel.reflectance;
            }

            if(!allow_transmission || fresnel.reflectance >= 1.0f) {
                return color(0.0f);
            }
            const auto transmitted = refract(sample_record.wi, vec3(0.0f, 0.0f, 1.0f),
                                             fresnel.cos_theta_t, eta_i, eta_t);
            if(!delta_aligned(transmitted, sample_record.wo)) {
                return color(0.0f);
            }
            auto value = specular_transmittance * (1.0f - fresnel.reflectance);
            if(ctx.mode == bsdf_context::transport_mode::RADIANCE) {
                const auto eta = eta_i / eta_t;
                value *= eta * eta;
            }
            return value;
        }

        global_precision pdf(const bsdf_context& ctx, const bsdf_sample_record& sample_record, measure meas) const override {
            if(meas != measure::discrete || sample_record.wi.z <= 0.0f) {
                return 0.0f;
            }

            const bool allow_reflection = (ctx.type_mask & reflection_lobe) != 0
                && (ctx.component == -1 || ctx.component == 0);
            const bool allow_transmission = (ctx.type_mask & transmission_lobe) != 0
                && (ctx.component == -1 || ctx.component == 1);

            const auto eta_i = sample_record.front_face ? ext_ior : int_ior;
            const auto eta_t = sample_record.front_face ? int_ior : ext_ior;
            const auto fresnel = fresnel_dielectric(sample_record.wi.z, eta_i, eta_t);

            if(sample_record.wi.z * sample_record.wo.z >= 0.0f) {
                if(!allow_reflection
                    || !delta_aligned(reflect(sample_record.wi, vec3(0.0f, 0.0f, 1.0f)), sample_record.wo)) {
                    return 0.0f;
                }
                return allow_transmission ? fresnel.reflectance : 1.0f;
            }

            if(!allow_transmission || fresnel.reflectance >= 1.0f) {
                return 0.0f;
            }
            const auto transmitted = refract(sample_record.wi, vec3(0.0f, 0.0f, 1.0f),
                                             fresnel.cos_theta_t, eta_i, eta_t);
            if(!delta_aligned(transmitted, sample_record.wo)) {
                return 0.0f;
            }
            return allow_reflection ? (1.0f - fresnel.reflectance) : 1.0f;
        }
    private:
        static constexpr unsigned int reflection_lobe = SPECULAR | REFLECTIVE | DIELECTRIC;
        static constexpr unsigned int transmission_lobe = SPECULAR | REFRACTIVE | DIELECTRIC;

        global_precision int_ior;
        global_precision ext_ior;
        color specular_reflectance;
        color specular_transmittance;
    };
}