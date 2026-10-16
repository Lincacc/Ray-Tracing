#pragma once

#include "integrator/photonmapper/photon_map.h"

#include <cstdint>
#include <vector>

namespace LinRT {
    // 一个工作单元的结果。clear() 只清长度，容量留给同一线程的下一段。
    class photon_batch {
    public:
        void clear() {
            photons.clear();
            boundaries.clear();
        }

        // 粒子开始时调用一次。这一段全部追踪完，再调用一次，用来封住最后一颗粒子。
        void next_particle() {
            boundaries.push_back(static_cast<std::uint32_t>(photons.size()));
        }

        void put(const photon& value) { photons.push_back(value); }

        std::vector<photon>& photon_buffer() { return photons; }

        std::size_t photon_count() const { return photons.size(); }

        std::size_t particle_count() const {
            return boundaries.empty() ? 0 : boundaries.size() - 1;
        }

        std::uint32_t particle_begin(std::size_t index) const { return boundaries[index]; }

        std::uint32_t particle_end(std::size_t index) const { return boundaries[index + 1]; }

        const photon& operator[](std::size_t index) const { return photons[index]; }

    private:
        std::vector<photon> photons;
        std::vector<std::uint32_t> boundaries;
    };
}