#pragma once

#include "emitter.h"
#include "geometry/hittable.h"

namespace LinRT {
    class area_emitter : public emitter {
    public:
        area_emitter(std::shared_ptr<hittable> shape, const color& radiance)
        : shape_ptr(std::move(shape)), radiance(radiance) {
            if(this->shape_ptr) {
                this->shape_ptr->set_emitter(this);
            }
        }

        bool sample(const vec3& shading_point, const vec2& sample, direct_sample& out) const override {
            if(!shape_ptr || shape_ptr->get_area() <= 0.0f) {
                return false;
            }

            vec3 p, n;
            global_precision pdf_area = 0.0f;
            shape_ptr->sample_position(sample, p, n, pdf_area);

            const auto to_light = p - shading_point;
            const auto dist2 = to_light.length_squared();
            if(dist2 < 1e-8f) {
                return false;
            }

            const auto dist = std::sqrt(dist2);
            const auto direction = to_light / dist;
            // 只朝法线一侧发光。着色点在背面时，这次抽样作废。
            const auto cos_light = dot(n, -direction);
            if(cos_light <= 0.0f) {
                return false;
            }

            out.point = p;
            out.normal = n;
            out.direction = direction;
            out.distance = dist;
            // 从面积密度换到立体角密度。积分器只除这一个 pdf，不要再乘 cos / 距离²。
            out.pdf = pdf_area * dist2 / cos_light;
            out.radiance = radiance;

            return true;
        }

        global_precision pdf(const vec3& shading_point, const vec3& direction) const override {
            if(!shape_ptr || shape_ptr->get_area() <= 0.0f) {
                return 0.0f;
            }

            const auto dir = normalize(direction);
            ray r(shading_point, dir);
            hit_record rec;
            // 和路径的 tmin 一致：积分器能打中的矩形，这里再问 pdf 也必须打中。
            if(!shape_ptr->hit(r, 0.001f, INFINITY, rec) || !rec.front_face) {
                return 0.0f;
            }

            // 正面命中时，rec.normal 就是该点朝外的法线。球体不能用形状上的一条公共法线。
            const auto cos_light = dot(rec.normal, -dir);
            if(cos_light <= 0.0f) {
                return 0.0f;
            }

            // dir 已归一化，rec.t 就是距离。
            return rec.t * rec.t / (shape_ptr->get_area() * cos_light);
        }

        // 飞出场景时看不见矩形，亮度只在命中发光面时由 emitted 给出。
        color eval(const vec3& direction) const override {
            return color(0.0f);
        }

        color emitted(const vec3& outward_normal, const vec3& outgoing) const override {
            if(dot(outward_normal, outgoing) <= 0.0f) {
                return color(0.0f);
            }

            return radiance;
        }

        bool is_infinite() const override {
            return false;
        }

    private:
        std::shared_ptr<hittable> shape_ptr;
        color radiance;
    };
}