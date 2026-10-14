#pragma once

#include "core/sobol.h"
#include "core/vec2.h"

#include <cstdint>

namespace LinRT {
    class sampler {
    public:
        // pixel_index 只决定扰乱，不改变样本编号。同一个像素的 128 个样本仍然是一条 Sobol 序列。
        void start(uint32_t pixel_index, uint32_t sample_index) {
            sample_index_ = sample_index;
            dimension_ = 0;
            pixel_seed_ = sobol::mix_bits(pixel_index);
        }

        global_precision next_1d() {
            const auto dimension = dimension_++;
            // 按像素、按维度各异或一次。相邻像素的灯上采样点不再是同一组圆。
            const auto scrambled = sobol::sample_bits(sample_index_, dimension) ^ scramble(dimension);
            return static_cast<global_precision>(scrambled) * (1.0f / 4294967296.0f);
        }

        vec2 next_2d() {
            return vec2(next_1d(), next_1d());
        }

    private:
        uint32_t scramble(uint32_t dimension) const {
            return sobol::mix_bits(pixel_seed_ + dimension * 0x9e3779b9u);
        }

        uint32_t sample_index_ = 0;
        uint32_t dimension_ = 0;
        uint32_t pixel_seed_ = 0;
    };
}