#pragma once

#include "integrator/photonmapper/photon_map.h"
#include "integrator/photonmapper/photon_gather.h"
#include "core/scene.h"
#include "core/sampler.h"
#include "core/onb.h"

#include <algorithm>
#include <cmath>

namespace LinRT {
    class photon_integrator {
    public:
        void preprocess(const scene& world) {
            const auto radius = world.bounding_sphere_radius();
            global_radius = 0.05f * radius;
            caustic_radius = 0.0125f * radius;

            auto caustic_requested = caustic_photons;
            if(!world.has_delta_bsdf()) {
                caustic_requested = 0;
            }

            photon_gather gather;
            const auto global = gather.run(
                world, photon_tracer::gather::surface, global_photons, max_depth - 1, rr_depth);
            if(global.full) {
                global_map.build(global.photons, 1.0f / static_cast<global_precision>(global.shot));
            }

            const auto caustic = gather.run(
                world, photon_tracer::gather::caustic, caustic_requested, max_depth - 1, rr_depth);
            if(caustic.full) {
                caustic_map.build(caustic.photons, 1.0f / static_cast<global_precision>(caustic.shot));
            }
        }

        color radiance(const ray& camera_ray, const scene& world, sampler& sampler) const {
            eye_query query;
            query.depth = 1;
            return li(camera_ray, world, sampler, query);
        }

    private:
        struct eye_query {
            int depth = 1;
            bool evaluated = false;
            bool missed = false;
            hit_record hit;
            bool emission = true;
            bool direct = true;
            bool indirect = true;
            bool caustic = true;
        };

        static constexpr int max_depth = 128;
        static constexpr int max_specular_depth = 4;
        static constexpr int rr_depth = 5;
        static constexpr int direct_samples = 16;
        static constexpr int glossy_samples = 32;
        static constexpr std::size_t global_photons = 250000;
        static constexpr std::size_t caustic_photons = 250000;
        static constexpr std::size_t lookup_size = 120;
        static constexpr unsigned int smooth_flags = DIFFUSE | GLOSSY;

        color li(const ray& camera_ray, const scene& world, sampler& sampler, eye_query query) const {
            if(!query.evaluated) {
                query.missed = !world.hit(camera_ray, 0.001f, INFINITY, query.hit);
                query.evaluated = true;
            }

            if(query.missed) {
                return query.emission ? world.eval_environment(camera_ray.direction) : color(0.0f);
            }

            color result(0.0f);
            const auto& rec = query.hit;
            if(query.emission && rec.emitter_ptr && rec.front_face) {
                result += rec.emitter_ptr->emitted(rec.normal, -normalize(camera_ray.direction));
            }

            if(query.depth >= max_depth || !rec.bsdf_ptr) {
                return result;
            }

            const auto flags = rec.bsdf_ptr->flags();
            const auto scattering = flags & (DIFFUSE | GLOSSY | SPECULAR);
            const auto diffuse = scattering == DIFFUSE;
            const auto specular = (flags & SPECULAR) != 0;
            const auto smooth = (flags & smooth_flags) != 0;
            const auto exhaustive = specular && query.depth < max_specular_depth;

            // 法线已经朝向入射一侧。背面不取光子。
            if(diffuse && dot(rec.normal, camera_ray.direction) < 0.0f) {
                const auto photon_depth = max_depth - query.depth;
                const auto reflectance = rec.bsdf_ptr->diffuse_reflectance() * (1.0f / M_PI<global_precision>);
                if(query.indirect && !global_map.empty()) {
                    result += global_map.estimate_irradiance(
                        rec.p, rec.normal, global_radius, photon_depth, lookup_size) * reflectance;
                }
                if(query.caustic && !caustic_map.empty()) {
                    result += caustic_map.estimate_irradiance(
                        rec.p, rec.normal, caustic_radius, photon_depth, lookup_size) * reflectance;
                }
            }

            if(exhaustive && query.indirect) {
                result += specular_lobes(camera_ray, rec, world, sampler, query.depth);
            }

            if(smooth && query.direct) {
                result += emitter_samples(camera_ray, rec, world, sampler, query.depth, diffuse, specular, exhaustive);
            }

            const auto sample_direct = query.direct && (smooth || (specular && !exhaustive));
            const auto sample_indirect = query.indirect && !diffuse && (smooth || (specular && !exhaustive));
            if(sample_direct || sample_indirect) {
                result += bsdf_samples(camera_ray, rec, world, sampler, query.depth, diffuse, specular, exhaustive,
                                       sample_direct, sample_indirect);
            }

            return result;
        }

        color specular_lobes(const ray& camera_ray, const hit_record& rec, const scene& world,
                             sampler& sampler, int depth) const {
            onb basis;
            basis.build_from_normal(rec.normal);
            color result(0.0f);

            for(int component = 0; component < rec.bsdf_ptr->component_count(); ++component) {
                if((rec.bsdf_ptr->component_flags(component) & SPECULAR) == 0) {
                    continue;
                }

                bsdf_context ctx;
                ctx.mode = bsdf_context::transport_mode::RADIANCE;
                ctx.component = component;
                bsdf_sample_record sample_record;
                prepare_sample(camera_ray, rec, basis, sample_record);

                // 分量已经指定，二维样本不参与选瓣。
                const auto weight = rec.bsdf_ptr->sample(ctx, sample_record, vec2(0.5f, 0.5f));
                if(weight.is_zero()) {
                    continue;
                }

                eye_query child;
                child.depth = depth + 1;
                result += weight * li(ray(rec.p, basis.local_to_world(sample_record.wo)), world, sampler, child);
            }

            return result;
        }

        color emitter_samples(const ray& camera_ray, const hit_record& rec, const scene& world, sampler& sampler,
                              int depth, bool diffuse, bool specular, bool exhaustive) const {
            const auto counts = sample_counts(depth, diffuse);
            onb basis;
            basis.build_from_normal(rec.normal);
            color result(0.0f);

            for(int i = 0; i < counts.emitter; ++i) {
                direct_sample light;
                if(!world.sample_emitter_direction(rec.p, sampler.next_2d(), light)) {
                    continue;
                }

                const auto shadow_origin = rec.p + rec.normal * 1e-4f;
                auto shadow_direction = light.direction;
                auto tmax = INFINITY;
                if(!std::isinf(light.distance)) {
                    const auto to_light = light.point - shadow_origin;
                    const auto distance = to_light.length();
                    if(distance < 1e-6f) {
                        continue;
                    }
                    shadow_direction = to_light / distance;
                    tmax = distance * (1.0f - 1e-4f);
                }

                bsdf_sample_record sample_record;
                prepare_sample(camera_ray, rec, basis, sample_record);
                sample_record.wo = basis.world_to_local(shadow_direction);

                bsdf_context ctx;
                ctx.mode = bsdf_context::transport_mode::RADIANCE;
                if(!specular || exhaustive) {
                    ctx.type_mask = smooth_flags;
                }
                const auto value = rec.bsdf_ptr->eval(ctx, sample_record, measure::solid_angle);
                if(value.is_zero()) {
                    continue;
                }

                hit_record occluder;
                const auto blocked = world.hit(ray(shadow_origin, shadow_direction), 0.001f, tmax, occluder)
                    && occluder.emitter_ptr != light.emitter_ptr;
                if(blocked) {
                    continue;
                }

                const auto light_pdf = world.pdf_emitter_direction(rec.p, shadow_direction);
                const auto bsdf_pdf = rec.bsdf_ptr->pdf(ctx, sample_record, measure::solid_angle);
                if(light_pdf <= 0.0f) {
                    continue;
                }

                const auto weight = power_heuristic(
                    light_pdf * static_cast<global_precision>(counts.emitter),
                    bsdf_pdf * static_cast<global_precision>(counts.bsdf)) * counts.emitter_weight;
                result += value * light.radiance * (weight / light_pdf);
            }

            return result;
        }

        color bsdf_samples(const ray& camera_ray, const hit_record& rec, const scene& world, sampler& sampler,
                           int depth, bool diffuse, bool specular, bool exhaustive,
                           bool sample_direct, bool sample_indirect) const {
            const auto counts = sample_counts(depth, diffuse);
            onb basis;
            basis.build_from_normal(rec.normal);
            color result(0.0f);

            for(int i = 0; i < counts.bsdf; ++i) {
                bsdf_context ctx;
                ctx.mode = bsdf_context::transport_mode::RADIANCE;
                ctx.lobe_sample = sampler.next_1d();
                if(!specular || exhaustive) {
                    ctx.type_mask = smooth_flags;
                }

                bsdf_sample_record sample_record;
                prepare_sample(camera_ray, rec, basis, sample_record);
                const auto weight = rec.bsdf_ptr->sample(ctx, sample_record, sampler.next_2d());
                if(weight.is_zero() || sample_record.pdf <= 0.0f) {
                    continue;
                }

                const ray scattered(rec.p, basis.local_to_world(sample_record.wo));
                hit_record occluder;
                const auto hit = world.hit(scattered, 0.001f, INFINITY, occluder);

                if(sample_direct) {
                    color emission(0.0f);
                    auto light_pdf = 0.0f;
                    if(hit && occluder.emitter_ptr && occluder.front_face) {
                        emission = occluder.emitter_ptr->emitted(occluder.normal, -normalize(scattered.direction));
                        light_pdf = world.pdf_emitter_direction(rec.p, scattered.direction);
                    } else if(!hit) {
                        emission = world.eval_environment(scattered.direction);
                        light_pdf = world.pdf_emitter_direction(rec.p, scattered.direction);
                    }

                    if(!emission.is_zero()) {
                        const auto mis = power_heuristic(
                            sample_record.pdf * static_cast<global_precision>(counts.bsdf),
                            light_pdf * static_cast<global_precision>(counts.emitter)) * counts.bsdf_weight;
                        result += emission * weight * mis;
                    }
                }

                if(sample_indirect) {
                    eye_query child;
                    child.depth = depth + 1;
                    child.evaluated = true;
                    child.missed = !hit;
                    child.hit = occluder;
                    child.emission = false;
                    result += weight * li(scattered, world, sampler, child) * counts.bsdf_weight;
                }
            }

            return result;
        }

        struct counts {
            int emitter = 1;
            int bsdf = 1;
            global_precision emitter_weight = 1.0f;
            global_precision bsdf_weight = 1.0f;
        };

        // 第一击才用 16 或 32。递归结果不会被像素直接看到，Mitsuba 在 depth > 1 时改成各 1 个样本。
        static counts sample_counts(int depth, bool diffuse) {
            counts value;
            if(depth == 1) {
                value.emitter = direct_samples;
                value.bsdf = diffuse ? direct_samples : glossy_samples;
                value.emitter_weight = 1.0f / static_cast<global_precision>(value.emitter);
                value.bsdf_weight = 1.0f / static_cast<global_precision>(value.bsdf);
            }
            return value;
        }

        static void prepare_sample(const ray& camera_ray, const hit_record& rec, const onb& basis,
                                   bsdf_sample_record& sample_record) {
            sample_record.p = rec.p;
            sample_record.n = normalize(rec.normal);
            sample_record.wi = basis.world_to_local(-normalize(camera_ray.direction));
            sample_record.front_face = rec.front_face;
        }

        static global_precision power_heuristic(global_precision pdf_strategy, global_precision pdf_other) {
            const auto a = pdf_strategy * pdf_strategy;
            const auto b = pdf_other * pdf_other;
            const auto sum = a + b;
            return sum > 0.0f ? a / sum : 0.0f;
        }

        photon_map global_map;
        photon_map caustic_map;
        global_precision global_radius = 0.0f;
        global_precision caustic_radius = 0.0f;
    };
}