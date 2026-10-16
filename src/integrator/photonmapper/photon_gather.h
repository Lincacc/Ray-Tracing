#pragma once

#include "accelerator/thread_pool.h"
#include "integrator/photonmapper/photon_batch.h"
#include "integrator/photonmapper/photon_tracer.h"

#include <algorithm>
#include <cstdint>
#include <mutex>
#include <thread>

namespace LinRT {
    class photon_gather {
    public:
        photon_tracer::result run(const scene& world, photon_tracer::gather kind, std::size_t requested_photons,
                                  int max_depth, int rr_depth) {
            photon_tracer::result outcome;
            if(requested_photons == 0 || max_depth < 1) {
                return outcome;
            }

            const auto workers = std::max(1u, std::thread::hardware_concurrency());
            // ParticleProcess：granularity 为 0 时，段长 = 目标光子数 / (16 × 线程数)。
            granularity = std::max<std::size_t>(1, requested_photons / (16 * static_cast<std::size_t>(workers)));
            requested = requested_photons;
            shot = 0;
            returned = 0;
            next_index = 0;
            photons.clear();
            photons.reserve(requested);

            {
                thread_pool pool(workers);
                for(unsigned worker = 0; worker < workers; ++worker) {
                    pool.submit([this, &world, kind, max_depth, rr_depth]{
                        // 每个线程持有一份 PhotonVector，clear() 后继续用同一块容量。
                        photon_batch batch;
                        sampler particle_sampler;
                        std::uint32_t begin = 0;
                        std::uint32_t end = 0;
                        while(generate_work(begin, end)) {
                            batch.clear();
                            for(std::uint32_t index = begin; index < end; ++index) {
                                batch.next_particle();
                                particle_sampler.start(0, index);
                                photon_tracer::trace_particle(
                                    world, kind, particle_sampler, max_depth, rr_depth, batch.photon_buffer());
                            }
                            batch.next_particle();
                            process_result(batch);
                        }
                    });
                }
                pool.wait();
            }

            outcome.photons = std::move(photons);
            outcome.shot = shot;
            outcome.full = outcome.photons.size() == requested;
            return outcome;
        }

    private:
        // 已返回的光子数达到目标，或发射超过 10 万仍几乎存不进，就不再发新段。
        bool generate_work(std::uint32_t& begin, std::uint32_t& end) {
            std::lock_guard<std::mutex> lock(mutex);
            if(returned >= requested
               || (shot > 100000 && unsuccessful(requested, photons.size(), shot))) {
                return false;
            }
            begin = next_index;
            end = next_index + static_cast<std::uint32_t>(granularity);
            next_index = end;
            return true;
        }

        // 按粒子顺序追加。某一颗放不进时，它已经计入 shot，同一段后面的粒子不计。
        void process_result(const photon_batch& batch) {
            std::lock_guard<std::mutex> lock(mutex);
            std::size_t counted = 0;
            for(std::size_t i = 0; i < batch.particle_count(); ++i) {
                ++counted;
                auto full = false;
                for(auto j = batch.particle_begin(i); j < batch.particle_end(i); ++j) {
                    if(photons.size() >= requested) {
                        full = true;
                        break;
                    }
                    photons.push_back(batch[static_cast<std::size_t>(j)]);
                }
                if(full) {
                    break;
                }
            }
            shot += counted;
            returned += batch.photon_count();
        }

        static bool unsuccessful(std::size_t needed, std::size_t stored, std::size_t shot_count) {
            return stored < needed && (stored == 0 || stored < shot_count / 1024);
        }

        std::mutex mutex;
        std::vector<photon> photons;
        std::size_t requested = 0;
        std::size_t granularity = 1;
        std::size_t shot = 0;
        std::size_t returned = 0;
        std::uint32_t next_index = 0;
    };
}