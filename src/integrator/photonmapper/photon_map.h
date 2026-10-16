#pragma once

#include "core/color.h"
#include "core/vec3.h"
#include "accelerator/aabb.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <queue>
#include <utility>
#include <vector>

namespace LinRT {
    // 方向是光子到达该点时的行进方向。Mitsuba 存的是 -toWorld(wi)，wi 指向来处，所以行进方向就是 ray.d。
    struct photon {
        vec3 position;
        vec3 direction;
        vec3 normal;
        color power;
        std::uint16_t depth = 0;
    };

    class photon_map {
    public:
        void build(std::vector<photon> gathered, global_precision scale) {
            photons = std::move(gathered);
            this->scale = scale;
            nodes.clear();
            if(!photons.empty()) {
                build_range(0, static_cast<int>(photons.size()));
            }
        }

        bool empty() const { return photons.empty(); }

        // max_depth 是 Li 传入的 m_maxDepth - depth。更深的光子不计入这次查询。
        color estimate_irradiance(const vec3& point, const vec3& normal, global_precision radius,
                                  int max_depth, std::size_t max_photons) const {
            if(photons.empty() || radius <= 0.0f || max_photons == 0) {
                return color(0.0f);
            }

            const auto squared_radius = radius * radius;
            const auto neighbors = nearest(point, squared_radius, max_photons);
            const auto inverse_squared_radius = 1.0f / squared_radius;
            color sum(0.0f);

            for(const auto& neighbor : neighbors) {
                const auto& candidate = photons[neighbor.second];
                if(static_cast<int>(candidate.depth) > max_depth) {
                    continue;
                }

                // wi 从交点指向光子来的方向。只收查询法线正面、且与光子法线同侧的贡献。
                const auto incoming = -candidate.direction;
                const auto incoming_dot_photon = dot(candidate.normal, incoming);
                const auto incoming_dot_shading = dot(normal, incoming);
                if(incoming_dot_shading <= 0.0f || dot(candidate.normal, normal) <= 1e-1f
                   || incoming_dot_photon <= 1e-2f) {
                    continue;
                }

                const auto square = 1.0f - neighbor.first * inverse_squared_radius;
                const auto weight = (incoming_dot_shading / incoming_dot_photon) * (square * square);
                sum += candidate.power * std::abs(weight);
            }

            return sum * (scale * 3.0f * (1.0f / M_PI<global_precision>) * inverse_squared_radius);
        }

    private:
        struct node {
            int left = -1;
            int right = -1;
            int start = 0;
            int count = 0;
            int axis = 0;
            global_precision split = 0.0f;
            bool leaf = false;
        };

        int build_range(int start, int count) {
            const auto index = static_cast<int>(nodes.size());
            nodes.push_back(node{});
            if(count <= 8) {
                nodes[static_cast<size_t>(index)].leaf = true;
                nodes[static_cast<size_t>(index)].start = start;
                nodes[static_cast<size_t>(index)].count = count;
                return index;
            }

            aabb box;
            for(int i = 0; i < count; ++i) {
                box.extend(photons[static_cast<size_t>(start + i)].position);
            }
            const auto extent = box.maximum - box.minimum;
            auto axis = 0;
            if(extent.y > extent[axis]) { axis = 1; }
            if(extent.z > extent[axis]) { axis = 2; }

            const auto mid = start + count / 2;
            std::nth_element(
                photons.begin() + start, photons.begin() + mid, photons.begin() + start + count,
                [axis](const photon& a, const photon& b) { return a.position[axis] < b.position[axis]; });

            const auto left = build_range(start, mid - start);
            const auto right = build_range(mid, start + count - mid);
            nodes[static_cast<size_t>(index)].axis = axis;
            nodes[static_cast<size_t>(index)].split = photons[static_cast<size_t>(mid)].position[axis];
            nodes[static_cast<size_t>(index)].left = left;
            nodes[static_cast<size_t>(index)].right = right;
            return index;
        }

        // 堆满之后才把搜索球收到当前最远邻居。核函数使用的仍是调用时的 squared_radius。
        std::vector<std::pair<global_precision, std::uint32_t>> nearest(
            const vec3& point, global_precision squared_radius, std::size_t max_photons) const {
            using neighbor = std::pair<global_precision, std::uint32_t>;
            std::priority_queue<neighbor> farthest;
            search(0, point, squared_radius, max_photons, farthest);

            std::vector<neighbor> result;
            result.reserve(farthest.size());
            while(!farthest.empty()) {
                result.push_back(farthest.top());
                farthest.pop();
            }
            return result;
        }

        void search(int index, const vec3& point, global_precision squared_radius, std::size_t max_photons,
                    std::priority_queue<std::pair<global_precision, std::uint32_t>>& farthest) const {
            const auto& current = nodes[static_cast<size_t>(index)];
            auto limit = squared_radius;
            if(farthest.size() == max_photons) {
                limit = std::min(limit, farthest.top().first);
            }

            if(current.leaf) {
                for(int i = 0; i < current.count; ++i) {
                    const auto photon_index = static_cast<std::uint32_t>(current.start + i);
                    const auto distance = (photons[photon_index].position - point).length_squared();
                    if(distance > limit) {
                        continue;
                    }
                    if(farthest.size() == max_photons) {
                        farthest.pop();
                    }
                    farthest.emplace(distance, photon_index);
                    if(farthest.size() == max_photons) {
                        limit = std::min(squared_radius, farthest.top().first);
                    }
                }
                return;
            }

            const auto delta = point[current.axis] - current.split;
            const auto near = delta <= 0.0f ? current.left : current.right;
            const auto far = delta <= 0.0f ? current.right : current.left;
            search(near, point, squared_radius, max_photons, farthest);
            if(delta * delta <= limit) {
                search(far, point, squared_radius, max_photons, farthest);
            }
        }

        std::vector<photon> photons;
        std::vector<node> nodes;
        global_precision scale = 1.0f;
    };
}