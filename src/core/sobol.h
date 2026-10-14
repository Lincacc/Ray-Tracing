#pragma once

#include "core.h"

#include <cstdint>
#include <array>

namespace LinRT {
    class sobol {
    public:
        // 返回 [0, 2^32) 的整数。dimension 从 0 开始。
        static uint32_t sample_bits(uint32_t index, uint32_t dimension) {
            const auto& directions = direction_table();
            if(dimension >= directions.size()) {
                // 表只覆盖前 16 维。更深的弹射退回哈希，避免再次落回同一条序列。
                return mix_bits(index + dimension * 0x9e3779b9u);
            }

            const auto gray = index ^ (index >> 1);
            uint32_t value = 0;
            const auto& row = directions[dimension];
            for(int bit = 0; bit < 32; ++bit) {
                if((gray >> bit) & 1u) {
                    value ^= row[static_cast<size_t>(bit)];
                }
            }
            return value;
        }

        static global_precision sampler(uint32_t index, uint32_t dimension) {
            return static_cast<global_precision>(sample_bits(index, dimension)) * (1.0f / 4294967296.0f);
        }

        static uint32_t mix_bits(uint32_t x) {
            x ^= x >> 16;
            x *= 0x7feb352du;
            x ^= x >> 15;
            x *= 0x846ca68bu;
            x ^= x >> 16;
            return x;
        }

    private:
        struct primitive {
            uint8_t degree;
            uint32_t a;
            uint32_t m[6];
        };

        // 1-based 维度 2..16 的本原多项式次数、系数 a、初始 m。
        static const std::array<primitive, 15>& primitives() {
            static const std::array<primitive, 15> table{{
                {1, 0u,  {1, 0, 0, 0, 0, 0}},
                {2, 1u,  {1, 3, 0, 0, 0, 0}},
                {3, 1u,  {1, 3, 1, 0, 0, 0}},
                {3, 2u,  {1, 1, 1, 0, 0, 0}},
                {4, 1u,  {1, 1, 3, 3, 0, 0}},
                {4, 4u,  {1, 3, 5, 13, 0, 0}},
                {5, 2u,  {1, 1, 5, 5, 17, 0}},
                {5, 4u,  {1, 1, 5, 5, 5, 0}},
                {5, 7u,  {1, 1, 7, 11, 19, 0}},
                {5, 11u, {1, 1, 5, 1, 1, 0}},
                {5, 13u, {1, 1, 1, 3, 11, 0}},
                {5, 14u, {1, 3, 5, 5, 13, 0}},
                {6, 1u,  {1, 3, 3, 9, 7, 49}},
                {6, 13u, {1, 1, 1, 15, 21, 21}},
                {6, 16u, {1, 3, 1, 13, 27, 49}},
            }};
            return table;
        }

        // directions[dimension][bit] 是第 bit+1 个方向数。dimension 0 是第一维。
        static const std::array<std::array<uint32_t, 32>, 16>& direction_table() {
            static const auto table = [] {
                std::array<std::array<uint32_t, 32>, 16> directions{};
                for(int bit = 0; bit < 32; ++bit) {
                    directions[0][static_cast<size_t>(bit)] = 1u << (31 - bit);
                }

                const auto& rows = primitives();
                for(size_t row = 0; row < rows.size(); ++row) {
                    const auto degree = rows[row].degree;
                    uint32_t m[33] = {};
                    for(uint8_t i = 1; i <= degree; ++i) {
                        m[i] = rows[row].m[i - 1];
                    }
                    for(int i = static_cast<int>(degree) + 1; i <= 32; ++i) {
                        uint32_t acc = (m[i - degree] << degree) ^ m[i - degree];
                        for(int k = 1; k < static_cast<int>(degree); ++k) {
                            if((rows[row].a >> (degree - 1 - k)) & 1u) {
                                acc ^= m[i - k] << k;
                            }
                        }
                        if(i < 32) {
                            acc &= (1u << i) - 1u;
                        }
                        m[i] = acc;
                    }
                    for(int i = 1; i <= 32; ++i) {
                        directions[row + 1][static_cast<size_t>(i - 1)] = m[i] << (32 - i);
                    }
                }
                return directions;
            }();
            return table;
        }
    };
}