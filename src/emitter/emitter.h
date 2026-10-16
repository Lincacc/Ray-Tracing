#pragma once

#include "core/color.h"
#include "core/ray.h"
#include "core/vec2.h"

#include <algorithm>

namespace LinRT {
    class emitter;

    struct direct_sample {
        vec3 direction;
        vec3 point;
        vec3 normal;
        global_precision pdf = 0.0f;
        // 无限远光源保持这个默认值。矩形光在 sample 里写成到着色点的距离。
        global_precision distance = INFINITY;
        color radiance;

        const emitter* emitter_ptr = nullptr;
    };

    // 从光源表面发出的一条射线。weight 是辐射度权重，不含位置密度和方向密度。
    // 面光源上：weight * pdf_position * pdf_direction = Le * cosθ。
    // selection_pdf 由 scene 写入；直接调用光源的 sample_ray 时它保持 1。
    struct emitter_ray_sample {
        ray emitted;
        vec3 normal;
        color weight;
        global_precision pdf_position = 0.0f;
        global_precision pdf_direction = 0.0f;
        global_precision selection_pdf = 1.0f;
        const emitter* emitter_ptr = nullptr;
    };

    class emitter {
    public:
        virtual ~emitter() = default;
        /**
         * @brief 从着色点采样一条指向光源的方向。采样失败返回 false。
         * @param shading_point 着色点
         * @param sample 采样点
         * @param out 采样结果
         * @return 采样成功返回 true，否则返回 false
         */
        virtual bool sample(const vec3& shading_point, const vec2& sample, direct_sample& out) const = 0;

        /**
         * @brief 给定世界空间方向，返回 sample_direct 会给出的立体角密度。打不中光源时返回 0。
         * @param shading_point 着色点
         * @param direction 方向
         * @return PDF
         */
        virtual global_precision pdf(const vec3& shading_point, const vec3& direction) const = 0;

        virtual color eval(const vec3& direction) const = 0;

        // 射线打到发光表面上时用。outgoing 从交点指向来时的方向。背面返回黑色。
        virtual color emitted(const vec3&, const vec3&) const { return color(0.0f); }

        virtual bool is_infinite() const = 0;

        // 下面四项是 Mitsuba AreaLight 的发射测度：面积上的位置，加上法线半球上的方向。
        // 环境光没有这个测度，对应函数返回 0 或 false。下一事件仍走上面的 sample / pdf。
        virtual bool sample_ray(const vec2& spatial, const vec2& directional, emitter_ray_sample& out) const = 0;

        virtual global_precision pdf_position() const = 0;

        // direction 不必归一化。背面返回 0。密度是立体角密度。
        virtual global_precision pdf_direction(const vec3& normal, const vec3& direction) const = 0;

        // 面光源功率 Le * π * area。环境光返回黑色。
        virtual color power() const = 0;
    };
}