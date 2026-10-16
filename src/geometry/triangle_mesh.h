#pragma once

#include "hittable.h"
#include "accelerator/sah_kdtree.h"

#include <algorithm>
#include <cstdint>
#include <vector>

namespace LinRT {
    class triangle_mesh : public hittable {
    public:
        // indices 每三个一组。corner_normals 与 indices 等长，每个角点一条法线。
        triangle_mesh(std::vector<vec3> positions, std::vector<std::uint32_t> indices, std::vector<vec3> corner_normals)
            : positions(std::move(positions)), indices(std::move(indices)), corner_normals(std::move(corner_normals)) {
            const auto triangle_count = this->indices.size() / 3;
            area_cdf.reserve(triangle_count);
            area = 0.0f;
            for(size_t i = 0; i < triangle_count; ++i) {
                const auto& p0 = this->positions[this->indices[3 * i]];
                const auto& p1 = this->positions[this->indices[3 * i + 1]];
                const auto& p2 = this->positions[this->indices[3 * i + 2]];
                area += 0.5f * cross(p1 - p0, p2 - p0).length();
                area_cdf.push_back(area);
            }

            build_accelerator();
        }

        bool hit(const ray& r, global_precision t_min, global_precision t_max, hit_record& rec) const override {
            bool hit_any = false;
            auto closest = t_max;
            // 叶子回调里若命中，必须把 closest 写成更近的 t。树用它裁掉更远的体素。
            kdtree.intersect(r, t_min, closest, [&](std::uint32_t triangle, global_precision mint, global_precision maxt) {
                if(hit_triangle(triangle, r, mint, maxt, closest, rec)) {
                    hit_any = true;
                }
            });

            return hit_any;
        }

        bool hit_triangle(std::uint32_t triangle, const ray& r, global_precision t0, global_precision t1,
                        global_precision& closest, hit_record& rec) const {
            const auto base = static_cast<size_t>(triangle) * 3;
            const auto& p0 = positions[indices[base]];
            const auto& p1 = positions[indices[base + 1]];
            const auto& p2 = positions[indices[base + 2]];
            const auto edge1 = p1 - p0;
            const auto edge2 = p2 - p0;
            const auto geometric = cross(edge1, edge2);
            if(geometric.length_squared() <= 0.0f) {
            return false;
            }

            const auto pvec = cross(r.direction, edge2);
            const auto det = dot(edge1, pvec);
            if(std::abs(det) < 1e-8f) {
            return false;
            }

            const auto inv_det = 1.0f / det;
            const auto tvec = r.origin - p0;
            const auto u = dot(tvec, pvec) * inv_det;
            if(u < 0.0f || u > 1.0f) {
            return false;
            }

            const auto qvec = cross(tvec, edge1);
            const auto v = dot(r.direction, qvec) * inv_det;
            if(v < 0.0f || u + v > 1.0f) {
            return false;
            }

            const auto t = dot(edge2, qvec) * inv_det;
            // t1 已经是这段体素远端和当前最近交点里更小的那个。
            if(t < t0 || t > t1) {
            return false;
            }

            closest = t;
            rec.t = t;
            rec.p = r.at(t);
            rec.bsdf_ptr = bsdf_ptr;
            rec.emitter_ptr = emitter_ptr;
            rec.set_face_normal(r, geometric / geometric.length());

            const auto b0 = 1.0f - u - v;
            auto shading = b0 * corner_normals[base] + u * corner_normals[base + 1] + v * corner_normals[base + 2];
            if(shading.length_squared() > 0.0f) {
                shading = shading / shading.length();
                if(dot(r.direction, shading) > 0.0f) {
                    shading = -shading;
                }
                rec.normal = shading;
            }

            return true;
        }

        aabb bounds() const override {
            aabb box;
            for(const auto& position : positions) {
                box.extend(position);
            }
            return box;
        }
        
        void sample_position(const vec2& sample, vec3& p, vec3& n, global_precision& pdf) const override {
            pdf = 0.0f;
            if(area <= 0.0f || area_cdf.empty()) {
                return;
            }

            // sample.x 按面积选三角形，余数再变成这块三角形上的第一个重心坐标。
            const auto target = std::min(sample.x, 0.999999f) * area;
            const auto it = std::upper_bound(area_cdf.begin(), area_cdf.end(), target);
            auto tri = static_cast<size_t>(it - area_cdf.begin());
            if(tri >= area_cdf.size()) {
                tri = area_cdf.size() - 1;
            }
            const auto prev = tri == 0 ? 0.0f : area_cdf[tri - 1];
            const auto span = area_cdf[tri] - prev;
            const auto u0 = span > 0.0f ? (target - prev) / span : 0.0f;

            const auto su = std::sqrt(u0);
            const auto b0 = 1.0f - su;
            const auto b1 = sample.y * su;
            const auto b2 = 1.0f - b0 - b1;

            const auto i0 = indices[3 * tri];
            const auto i1 = indices[3 * tri + 1];
            const auto i2 = indices[3 * tri + 2];
            p = b0 * positions[i0] + b1 * positions[i1] + b2 * positions[i2];
            const auto geometric = cross(positions[i1] - positions[i0], positions[i2] - positions[i0]);
            n = geometric / geometric.length();
            pdf = 1.0f / area;
        }

        void translate(const vec3& delta) {
            for(auto& position : positions) {
                position += delta;
            }

            build_accelerator();
        }

        void build_accelerator() {
            const auto triangle_count = indices.size() / 3;
            std::vector<aabb> boxes(triangle_count);
            for(size_t i = 0; i < triangle_count; ++i) {
                boxes[i].extend(positions[indices[3 * i]]);
                boxes[i].extend(positions[indices[3 * i + 1]]);
                boxes[i].extend(positions[indices[3 * i + 2]]);
            }
            kdtree.build(boxes);
        }
    private:
        sah_kdtree kdtree;
        std::vector<vec3> positions;
        std::vector<std::uint32_t> indices;
        std::vector<vec3> corner_normals;
        std::vector<global_precision> area_cdf;
    };
}