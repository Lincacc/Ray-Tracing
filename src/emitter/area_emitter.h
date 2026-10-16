#pragma once

#include "emitter.h"
#include "geometry/hittable.h"
#include "core/onb.h"

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

        // Mitsuba AreaLight::sampleRay：余弦加权半球，返回 m_power = Le * π * area。
        // 余弦密度是 cosθ/π，所以 weight * pdf_position * pdf_direction = Le * cosθ。
        bool sample_ray(const vec2& spatial, const vec2& directional, emitter_ray_sample& out) const override {
            if(!shape_ptr || shape_ptr->get_area() <= 0.0f) {
                return false;
            }

            vec3 p, n;
            global_precision pdf_area = 0.0f;
            shape_ptr->sample_position(spatial, p, n, pdf_area);
            if(pdf_area <= 0.0f) {
                return false;
            }

            // 和 diffuse_bsdf 相同的余弦半球：cosθ = sqrt(ξ0)，φ = 2π ξ1。密度 cosθ/π。
            const auto cos_theta = std::sqrt(directional.x);
            const auto sin_theta = std::sqrt(1.0f - directional.x);
            const auto phi = 2.0f * M_PI<global_precision> * directional.y;
            const vec3 local(std::cos(phi) * sin_theta, std::sin(phi) * sin_theta, cos_theta);

            onb frame;
            frame.build_from_normal(n);

            out.emitted = ray(p, frame.local_to_world(local));
            out.normal = n;
            out.pdf_position = pdf_area;
            out.pdf_direction = cos_theta * (1.0f / M_PI<global_precision>);
            out.weight = power();
            out.selection_pdf = 1.0f;
            out.emitter_ptr = this;
            return true;
        }

        global_precision pdf_position() const override {
            if(!shape_ptr || shape_ptr->get_area() <= 0.0f) {
                return 0.0f;
            }
            return 1.0f / shape_ptr->get_area();
        }

        global_precision pdf_direction(const vec3& normal, const vec3& direction) const override {
            const auto n = normalize(normal);
            const auto d = normalize(direction);
            const auto cos_theta = dot(n, d);
            if(cos_theta <= 0.0f) {
                return 0.0f;
            }
            return cos_theta * (1.0f / M_PI<global_precision>);
        }

        color power() const override {
            if(!shape_ptr || shape_ptr->get_area() <= 0.0f) {
                return color(0.0f);
            }
            // Mitsuba AreaLight::setParent：m_power = m_radiance * π * surfaceArea。
            return radiance * (M_PI<global_precision> * shape_ptr->get_area());
        }

    private:
        std::shared_ptr<hittable> shape_ptr;
        color radiance;
    };
}