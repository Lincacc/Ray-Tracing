#pragma once

#include "hittable.h"

namespace LinRT {
    class rectangle : public hittable {
    public:
        // corner 是一个角点，edge_u、edge_v 是从它出发的两条边，参数范围都是 [0, 1]。
        rectangle(const vec3 corner, const vec3& edge_u, const vec3& edge_v) 
        : corner(corner), edge_u(edge_u), edge_v(edge_v) {
            const auto n = cross(edge_u, edge_v);
            area = n.length();
            if(area > 0.0f) {
                normal = n / area;
            }
        }

        bool hit(const ray& r, global_precision t_min, global_precision t_max, hit_record& rec) const override {
            const auto denom = dot(normal, r.direction);
            if(std::abs(denom) < 1e-8f) {
                return false;
            }

            const auto t = dot(corner - r.origin, normal) / denom;
            if(t < t_min || t > t_max) {
                return false;
            }

            const auto p = r.at(t);
            const auto d = p - corner;
            const auto uu = dot(edge_u, edge_u);
            const auto vv = dot(edge_v, edge_v);
            const auto uv = dot(edge_u, edge_v);
            const auto du = dot(d, edge_u);
            const auto dv = dot(d, edge_v);
            const auto det = uu * vv - uv * uv;
            if(std::abs(det) < 1e-12f) {
                return false;
            }

            const auto u = (vv * du - uv * dv) / det;
            const auto v = (uu * dv - uv * du) / det;
            if(u < 0.0f || u > 1.0f || v < 0.0f || v > 1.0f) {
                return false;
            }

            rec.t = t;
            rec.p = p;
            rec.normal = normal;
            rec.bsdf_ptr = nullptr;
            rec.emitter_ptr = emitter_ptr;
            rec.set_face_normal(r, normal);

            return true;
        }

        aabb bounds() const override {
            aabb box;
            box.extend(corner);
            box.extend(corner + edge_u);
            box.extend(corner + edge_v);
            box.extend(corner + edge_u + edge_v);
            
            return box;
        }

        void sample_position(const vec2& sample, vec3& p, vec3& n, global_precision& pdf) const override {
            p = corner + sample.x * edge_u + sample.y * edge_v;
            n = normal;
            pdf = 1.0f / area;
        }

    private:
        vec3 corner;
        vec3 edge_u;
        vec3 edge_v;
        vec3 normal;
    };
}