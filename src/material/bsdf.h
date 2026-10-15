#pragma once 

#include "core/ray.h"
#include "core/color.h"
#include "core/vec2.h"

namespace LinRT {
    enum bsdf_flags : unsigned int {
        NONE        = 0,
        DIFFUSE     = 1u << 0,
        SPECULAR    = 1u << 1,
        REFLECTIVE  = 1u << 2,
        TRANSPARENT = 1u << 3,
        GLOSSY      = 1u << 4,
        METALLIC    = 1u << 5,
        DIELECTRIC  = 1u << 6,
        REFRACTIVE  = 1u << 7
    };

    struct bsdf_context {
        enum class transport_mode {
            RADIANCE,
            IMPORTANCE,
        };

        transport_mode mode = transport_mode::RADIANCE;
        unsigned int type_mask = 0xffffffffu;
        int component = -1;
        global_precision lobe_sample = 0.0f;
    };

    struct bsdf_sample_record {
        vec3 wi; // 入射方向，局部坐标，从交点指向外
        vec3 wo; // 出射方向，局部坐标，从交点指向外
        vec3 p;
        vec3 n;
        vec3 uv;
        global_precision eta = 1.0f;
        bool front_face = true; // true：射线在外部介质。玻璃据此选择 eta_i、eta_t
        global_precision pdf = 0.0f;     // 本次 sample 抽到 wo 的立体角密度
        unsigned int sampled_type = NONE; // 本次抽中的分量类型
        int sampled_component = -1;
    };

    class bsdf {
        // 所有计算均在局部坐标系下进行        
    public:
        virtual unsigned int flags() const = 0;

        /**
         * @brief 1. 采样：根据入射方向和随机数，生成一条出射光线，并返回加权颜色 (f * cos_theta / pdf)
         * @param bRec 采样记录，传入 wi，输出 wo
         * @param sample 2D 随机数（用于支撑 QMC / Sobol 采样）
         * @return color 经过余弦修正和 PDF 归一化的权重值
         */
        virtual color sample(const bsdf_context& ctx, bsdf_sample_record& sample_record, const vec2& sample) const = 0;

        /**
         * @brief 2. 求值：给定入射方向 wi 和出射方向 wo，计算表面反射率 f * cos_theta
         */
        virtual color eval(const bsdf_context& ctx, const bsdf_sample_record& sample_record) const = 0;

        /**
         * @brief 3. 概率密度：给定 wi 和 wo，计算通过该材质采样到 wo 的概率密度函数值 (PDF)
         */
        virtual global_precision pdf(const bsdf_context& ctx, const bsdf_sample_record& sample_record) const = 0;
    };
}