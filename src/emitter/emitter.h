#pragma once

#include "core/color.h"
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

    class emitter {
    public:
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
    };
}