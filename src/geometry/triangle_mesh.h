#pragma once

#include "hittable.h"

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
        }

        bool hit(const ray& r, global_precision t_min, global_precision t_max, hit_record& rec) const override {
            bool hit_any = false;
            auto closest = t_max;

            for(size_t i = 0; i < indices.size(); i += 3) {
                const auto& p0 = positions[indices[i]];
                const auto& p1 = positions[indices[i + 1]];
                const auto& p2 = positions[indices[i + 2]];
                const auto edge1 = p1 - p0;
                const auto edge2 = p2 - p0;
                const auto geometric = cross(edge1, edge2);
                if(geometric.length_squared() <= 0.0f) {
                    continue;
                }

                // Möller–Trumbore。u、v 是 p1、p2 的重心坐标，p0 的坐标是 1-u-v。
                const auto pvec = cross(r.direction, edge2);
                const auto det = dot(edge1, pvec);
                if(std::abs(det) < 1e-8f) {
                    continue;
                }
                const auto inv_det = 1.0f / det;
                const auto tvec = r.origin - p0;
                const auto u = dot(tvec, pvec) * inv_det;
                if(u < 0.0f || u > 1.0f) {
                    continue;
                }
                const auto qvec = cross(tvec, edge1);
                const auto v = dot(r.direction, qvec) * inv_det;
                if(v < 0.0f || u + v > 1.0f) {
                    continue;
                }
                const auto t = dot(edge2, qvec) * inv_det;
                if(t < t_min || t > closest) {
                    continue;
                }

                closest = t;
                hit_any = true;
                rec.t = t;
                rec.p = r.at(t);
                rec.bsdf_ptr = bsdf_ptr;
                rec.emitter_ptr = emitter_ptr;
                rec.set_face_normal(r, geometric / geometric.length());

                const auto b0 = 1.0f - u - v;
                auto shading = b0 * corner_normals[i] + u * corner_normals[i + 1] + v * corner_normals[i + 2];
                if(shading.length_squared() > 0.0f) {
                    shading = shading / shading.length();
                    // 着色法线翻到射线进来的一侧，局部空间里 wi.z 才是正的。front_face 仍是几何法线的结果。
                    if(dot(r.direction, shading) > 0.0f) {
                        shading = -shading;
                    }
                    rec.normal = shading;
                }
            }

            return hit_any;
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
        }

    private:
        std::vector<vec3> positions;
        std::vector<std::uint32_t> indices;
        std::vector<vec3> corner_normals;
        std::vector<global_precision> area_cdf;
    };
}