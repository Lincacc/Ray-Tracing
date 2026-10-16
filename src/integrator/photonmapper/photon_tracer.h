#pragma once

#include "integrator/photonmapper/photon_map.h"
#include "core/scene.h"
#include "core/sampler.h"
#include "core/onb.h"

#include <algorithm>
#include <cstdint>
#include <vector>

namespace LinRT {
    class photon_tracer {
    public:
        enum class gather { surface, caustic };

        struct result {
            std::vector<photon> photons;
            std::size_t shot = 0;
            bool full = false;
        };

        result gather_map(const scene& world, gather kind, std::size_t requested, int max_depth, int rr_depth) const {
            result gathered;
            if(requested == 0 || max_depth < 1) {
                return gathered;
            }

            gathered.photons.reserve(requested);
            for(std::uint32_t index = 0; gathered.photons.size() < requested; ++index) {
                if(gathered.shot > 100000 && unsuccessful(requested, gathered.photons.size(), gathered.shot)) {
                    break;
                }

                sampler particle_sampler;
                particle_sampler.start(0, index);
                trace_particle(world, kind, particle_sampler, max_depth, rr_depth, requested, gathered);
                ++gathered.shot;
            }

            gathered.full = gathered.photons.size() == requested;
            return gathered;
        }

    private:
        static bool unsuccessful(std::size_t needed, std::size_t stored, std::size_t shot) {
            return stored < needed && (stored == 0 || stored < shot / 1024);
        }

        static void trace_particle(const scene& world, gather kind, sampler& particle_sampler,
                                   int max_depth, int rr_depth, std::size_t requested, result& gathered) {
            emitter_ray_sample emitted;
            if(!world.sample_emitter_ray(particle_sampler.next_2d(), particle_sampler.next_2d(), emitted)) {
                return;
            }

            ray current = emitted.emitted;
            color throughput(1.0f);
            auto delta = false;

            // depth 从 1 开始。第一击不存储，直接光留给眼睛路径的光源采样。
            for(int depth = 1; !throughput.is_zero() && depth <= max_depth; ) {
                hit_record rec;
                if(!world.hit(current, 0.001f, INFINITY, rec) || !rec.bsdf_ptr) {
                    break;
                }

                store_if_accepted(kind, rec, current, throughput * emitted.weight, depth, delta, requested, gathered);

                onb basis;
                basis.build_from_normal(rec.normal);
                bsdf_context ctx;
                ctx.mode = bsdf_context::transport_mode::IMPORTANCE;
                ctx.lobe_sample = particle_sampler.next_1d();

                bsdf_sample_record sample_record;
                sample_record.p = rec.p;
                sample_record.n = normalize(rec.normal);
                sample_record.wi = basis.world_to_local(-normalize(current.direction));
                sample_record.front_face = rec.front_face;

                const auto weight = rec.bsdf_ptr->sample(ctx, sample_record, particle_sampler.next_2d());
                if(weight.is_zero()) {
                    break;
                }

                const auto outgoing = basis.local_to_world(sample_record.wo);
                const auto incoming = basis.local_to_world(sample_record.wi);
                // 几何法线和着色法线半球不一致时丢弃。这里两者是同一条法线，全内反射的折射仍能通过。
                if(dot(rec.normal, incoming) * sample_record.wi.z <= 0.0f
                   || dot(rec.normal, outgoing) * sample_record.wo.z <= 0.0f) {
                    break;
                }

                throughput *= weight;
                delta = (sample_record.sampled_type & SPECULAR) != 0;
                current = ray(rec.p, outgoing);

                if(depth++ >= rr_depth) {
                    const auto survival = std::min(std::max({throughput.r, throughput.g, throughput.b}), 0.95f);
                    if(particle_sampler.next_1d() >= survival) {
                        break;
                    }
                    throughput /= survival;
                }
            }
        }

        // 上一跳是 delta 才进焦散图。光滑玻璃自己没有 DIFFUSE/GLOSSY，光子不存在玻璃上。
        static void store_if_accepted(gather kind, const hit_record& rec, const ray& incoming, const color& power,
                                      int depth, bool delta, std::size_t requested, result& gathered) {
            const auto flags = rec.bsdf_ptr->flags();
            const auto smooth_reflection = (flags & DIFFUSE) != 0 || (flags & GLOSSY) != 0;
            const auto accepted = smooth_reflection && depth > 1 && (kind == gather::caustic ? delta : !delta);
            if(!accepted || gathered.photons.size() >= requested) {
                return;
            }

            photon stored;
            stored.position = rec.p;
            stored.direction = normalize(incoming.direction);
            stored.normal = rec.normal;
            stored.power = power;
            stored.depth = static_cast<std::uint16_t>(depth);
            gathered.photons.push_back(stored);
        }
    };
}