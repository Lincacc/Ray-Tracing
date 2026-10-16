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

        static void trace_particle(const scene& world, gather kind, sampler& particle_sampler,
                                int max_depth, int rr_depth, std::vector<photon>& out) {
            emitter_ray_sample emitted;
            if(!world.sample_emitter_ray(particle_sampler.next_2d(), particle_sampler.next_2d(), emitted)) {
                return;
            }

            ray current = emitted.emitted;
            color throughput(1.0f);
            auto delta = false;

            for(int depth = 1; !throughput.is_zero() && depth <= max_depth; ) {
                hit_record rec;
                if(!world.hit(current, 0.001f, INFINITY, rec) || !rec.bsdf_ptr) {
                    break;
                }

                store_if_accepted(kind, rec, current, throughput * emitted.weight, depth, delta, out);

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

        static void store_if_accepted(gather kind, const hit_record& rec, const ray& incoming, const color& power,
                                    int depth, bool delta, std::vector<photon>& out) {
            const auto flags = rec.bsdf_ptr->flags();
            const auto smooth_reflection = (flags & DIFFUSE) != 0 || (flags & GLOSSY) != 0;
            const auto accepted = smooth_reflection && depth > 1 && (kind == gather::caustic ? delta : !delta);
            if(!accepted) {
                return;
            }

            photon stored;
            stored.position = rec.p;
            stored.direction = normalize(incoming.direction);
            stored.normal = rec.normal;
            stored.power = power;
            stored.depth = static_cast<std::uint16_t>(depth);
            out.push_back(stored);
        }        

    private:
        static bool unsuccessful(std::size_t needed, std::size_t stored, std::size_t shot) {
            return stored < needed && (stored == 0 || stored < shot / 1024);
        }
    };
}