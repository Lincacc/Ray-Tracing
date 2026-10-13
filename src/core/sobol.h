#pragma once

#include "core.h"

#include <cstdint>

namespace LinRT {
    class sobol {
    public:
        static uint32_t bit_reverse(uint32_t n) {
            n = (n << 16) | (n >> 16);
            n = ((n & 0x00ff00ff) << 8) | ((n & 0xff00ff00) >> 8);
            n = ((n & 0x0f0f0f0f) << 4) | ((n & 0xf0f0f0f0) >> 4);
            n = ((n & 0x33333333) << 2) | ((n & 0xcccccccc) >> 2);
            n = ((n & 0x55555555) << 1) | ((n & 0xaaaaaaaa) >> 1);
            return n;
        }

        /**
         * @brief 生成 Sobol 序列采样点
         * 
         * @param index     当前采样点的索引 (从 0 开始)
         * @param dimension 当前采样点的维度 (从 1 开始)
         * @return global_precision 采样点的值 (范围: [0, 1))
         */
        static global_precision sampler(uint32_t index, uint32_t dimension) {
            // 一维：标准基数 2 的 Van der Corput 序列
            if(dimension == 1) {
                return static_cast<global_precision>(bit_reverse(index)) * (1.0f / 4294967296.0f);
            }

            // 多维度：基于 Gray 码与方向数的计算
            uint32_t v = 0;
            uint32_t gray = index ^ (index >> 1);

            for(int i=0; i<32; ++i) {
                if((gray >> i) & 1) {
                    v ^= get_direction_number(dimension, i);
                }
            }

            return static_cast<global_precision>(v) * (1.0f / 4294967296.0f);
        }

    private:
        // 获取指定维度的预定义方向数
        static uint32_t get_direction_number(uint32_t dimension, uint32_t index) {
            static const uint32_t dir_dim1[8] = {
                0x80000000, 0x40000000, 0x20000000, 0x10000000,
                0x08000000, 0x04000000, 0x02000000, 0x01000000
            };

            static const uint32_t dir_dim2[8] = {
                0x80000000, 0xc0000000, 0xa0000000, 0x30000000,
                0x68000000, 0x04000000, 0x52000000, 0x0b000000
            };

            if(dimension == 2) {
                return (index < 8) ? dir_dim1[index] : (1u << (31 - index));
            } else if(dimension == 3) {
                return (index < 8) ? dir_dim2[index] : (1u << (31 - index));
            }
            
            return (1u << (31 - index));
        }
    };
}