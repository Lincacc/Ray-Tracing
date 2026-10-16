#pragma once

#include "core/ray.h"
#include "material/bsdf.h"
#include "emitter/emitter.h"

#include <memory>

namespace LinRT {
    class bsdf;

    struct hit_record {
        global_precision t;
        vec3 p;
        vec3 normal;
        bool front_face;
        std::shared_ptr<bsdf> bsdf_ptr;
        const emitter* emitter_ptr = nullptr;

        // 统一法线方向：确保法线永远指向光线入射方向的反方向
        inline void set_face_normal(const ray& r, const vec3& outward_normal) {
            front_face = dot(r.direction, outward_normal) < 0;
            normal = front_face ? outward_normal : -outward_normal;
        }
    };

    class hittable {
    public:
        /**
         * @brief 检测光线与物体的交点
         * 
         * @param r 光线
         * @param t_min 光线参数的最小值
         * @param t_max 光线参数的最大值
         * @param rec 交点信息记录
         */
        virtual bool hit(const ray& r, global_precision t_min, global_precision t_max, hit_record& rec) const = 0;
    
        /**
         * @brief 设置物体的材质
         * 
         * @param mat 材质
         */
        virtual void set_material(std::shared_ptr<bsdf> bsdf) { bsdf_ptr = bsdf; }

        /**
         * @brief 获取物体的材质
         * 
         * @return 材质
         */
        virtual std::shared_ptr<bsdf> get_material() const { return bsdf_ptr; }

        virtual global_precision get_area() const { return area; }

        virtual void set_emitter(const emitter* emitter) { emitter_ptr = emitter; }

        virtual void sample_position(const vec2& sample, vec3& p, vec3& n, global_precision& pdf) const = 0;

    protected:
        std::shared_ptr<bsdf> bsdf_ptr;
        global_precision area = 0.0f;
        const emitter* emitter_ptr = nullptr;
    };
}