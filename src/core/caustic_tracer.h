#pragma once

#include "core/camera.h"
#include "core/hittable_list.h"
#include "core/sampler.h"
#include "emitter/emitter.h"
#include "material/bsdf.h"

#include <vector>

namespace LinRT {
    inline void splat_caustic(const hittable_list& world,
                            const std::vector<std::shared_ptr<emitter>>& emitters,
                            const camera& camera,
                            std::vector<color>& pixels,
                            int width, int height, int samples_per_pix) {
        const auto light_count = static_cast<int>(emitters.size());
        if(light_count <= 0 || width <= 1 || height <= 1) {
            return;
        }

        const auto paths = static_cast<long long>(width) * height * samples_per_pix;
        const auto film = camera.film_area();
        if(film <= 0.0f) {
            return;
        }

        for(long long i = 0; i < paths; ++i) {
            sampler sampler;
            sampler.start(static_cast<uint32_t>(i), 0xA341316Cu);

            const auto pick = sampler.next_1d() * static_cast<global_precision>(light_count);
            auto index = static_cast<int>(pick);
            if(index >= light_count) {
                index = light_count - 1;
            }

            emitter_ray_sample emitted;
            if(!emitters[static_cast<size_t>(index)]->sample_ray(sampler.next_2d(), sampler.next_2d(), emitted)) {
                continue;
            }

            // 均匀选灯的概率是 1/n，beta 里要除掉它。

            auto beta = emitted.beta * static_cast<global_precision>(light_count);
            auto bounce = emitted.emitted;
            bool specular_bounce = false;

            for(int depth = 0; depth < 8; ++depth) {
                hit_record rec;
                if(!world.hit(bounce, 0.001f, INFINITY, rec) || !rec.bsdf_ptr) {
                    break;
                }

                const auto flags = rec.bsdf_ptr->flags();
                const auto delta = (flags & SPECULAR) !=0 && (flags & (DIFFUSE | GLOSSY)) == 0;
                if(delta) {
                    onb basis;
                    basis.build_from_normal(rec.normal);

                    bsdf_sample_record sample_record;
                    sample_record.wi = basis.world_to_local(-normalize(bounce.direction));
                    sample_record.front_face = rec.front_face;

                    bsdf_context ctx;
                    ctx.mode = bsdf_context::transport_mode::RADIANCE;
                    const auto weight = rec.bsdf_ptr->sample(ctx, sample_record, sampler.next_2d());
                    if(weight.is_zero()) {
                        break;
                    }

                    // 球体一进一出，两次 (η_i/η_t)² 乘起来是 1，和眼睛路径相同。
                    beta *= weight;
                    bounce = ray(rec.p, basis.local_to_world(sample_record.wo));
                    specular_bounce = true;
                    continue;
                }

                if(!specular_bounce || (flags & DIFFUSE) == 0) {
                    break;
                }

                global_precision s = 0.0f, t = 0.0f, cos_camera = 0.0f;
                if(!camera.project(rec.p, s, t, cos_camera)) {
                    break;
                }

                const auto to_camera = camera.origin - rec.p;
                const auto dist2 = to_camera.length_squared();
                const auto dist = std::sqrt(dist2);
                if(dist <= 1e-4f) {
                    break;
                }

                onb basis;
                basis.build_from_normal(rec.normal);
                
                bsdf_sample_record sample_record;
                sample_record.wi = basis.world_to_local(-normalize(bounce.direction));
                sample_record.wo = basis.world_to_local(to_camera / dist);
                if(sample_record.wo.z <= 0.0f) {
                    break;
                }

                bsdf_context ctx;
                const auto f_cos = rec.bsdf_ptr->eval(ctx, sample_record);
                if(f_cos.is_zero()) {
                    break;
                }

                hit_record blocker;
                const auto shadow_dir = to_camera / dist;
                const ray shadow_ray(rec.p + rec.normal * 1e-4f, shadow_dir);
                if(world.hit(shadow_ray, 0.001f, dist * (1.0f - 1e-4f), blocker)) {
                    break;
                }

                // pbrt ConnectToCamera：beta 已除掉灯的密度，再除掉该像素的立体角。
                const auto contrib = beta * f_cos
                    / (film * cos_camera * cos_camera * cos_camera * dist2);
                const auto x = std::min(width - 1, static_cast<int>(s * static_cast<global_precision>(width - 1)));
                const auto y = std::min(height - 1, static_cast<int>(t * static_cast<global_precision>(height - 1)));
                pixels[static_cast<size_t>(y) * width + x] += contrib / static_cast<global_precision>(samples_per_pix);
                break;
            }
        }
    }
}