#pragma once 

#include "hittable.h"

namespace LinRT {
    class sphere : public hittable {
    public:
        sphere() = default;
        sphere(const vec3& center, global_precision radius) : center(center), radius(radius) {
            area = 4.0f * M_PI<global_precision> * radius * radius;
        }

        /**
         * @brief 检测光线与球体的交点
         * 
         * @param r 光线
         * @param t_min 光线参数的最小值
         * @param t_max 光线参数的最大值
         * @param rec 交点信息记录
         */
        virtual bool hit(const ray& r, global_precision t_min, global_precision t_max, hit_record& rec) const override {
            vec3 oc = r.origin - center;
            auto a = r.direction.length_squared();
            auto half_b = dot(oc, r.direction);
            auto c = oc.length_squared() - radius * radius;

            auto discriminant = half_b * half_b - a * c;
            if (discriminant < 0) return false;
            auto sqrtd = std::sqrt(discriminant);

            // 寻找最近的有效正根
            auto root = (-half_b - sqrtd) / a;
            if (root < t_min || root > t_max) {
                root = (-half_b + sqrtd) / a;
                if (root < t_min || root > t_max)
                    return false;
            }

            rec.t = root;
            rec.p = r.at(rec.t);
            rec.bsdf_ptr = bsdf_ptr;
            rec.emitter_ptr = emitter_ptr;

            // 计算向外的法线并统一朝向
            vec3 outward_normal = (rec.p - center) / radius;
            rec.set_face_normal(r, outward_normal);

            return true;
        }

        void sample_position(const vec2& sample, vec3& p, vec3& n, global_precision& pdf) const override {
            const auto z = 1.0f - 2.0f * sample.y;
            const auto ring = std::sqrt(std::max(0.0f, 1.0f - z * z));
            const auto phi = 2.0f * M_PI<global_precision> * sample.x;

            n = vec3(ring * std::cos(phi), ring * std::sin(phi), z);
            p = center + radius * n;
            pdf = 1.0f / area;
        }

    public:
        vec3 center;
        global_precision radius;
    };
}