#pragma once

#include "core/hittable_list.h"
#include "core/color.h"
#include "core/onb.h"
#include "emitter/env_emitter.h"
#include "core/sobol.h"

namespace LinRT {
    constexpr int MAX_DEPTH = 50;

    color trace(ray& r, const hittable_list& world, std::shared_ptr<emitter> emitter, uint32_t sample_index) {
        constexpr int rr_depth = 5;
        global_precision eta = 1.0f;
        color radiance(0.0f);
        color throughput(1.0f);
        for(int depth = 0; depth < MAX_DEPTH; ++depth) {
            hit_record rec;

            if(!world.hit(r, 0.001f, INFINITY, rec)) {
                radiance += throughput * emitter->eval(r.direction);
                break;
            }

            if(!rec.bsdf_ptr) {
                break;
            }

            onb onb;
            onb.build_from_normal(rec.normal);

            bsdf_sample_record sample_record;
            sample_record.p = rec.p;
            sample_record.n = normalize(rec.normal);
            // mitsuba: wi 指向上一跳的来源方向
            sample_record.wi = onb.world_to_local(-normalize(r.direction));

            // 随机 2D 采样数
            // 维度 0、1 留给像素。每次弹射占用 3 个维度：BSDF 两个，俄罗斯轮盘一个。
            auto sample_2d = vec2(
                sobol::sampler(sample_index, 2 + depth * 3), 
                sobol::sampler(sample_index, 3 + depth * 3));

            // 采样 bsdf
            bsdf_context ctx;
            ctx.mode = bsdf_context::transport_mode::RADIANCE;
            auto bsdf_weight = rec.bsdf_ptr->sample(ctx, sample_record, sample_2d);
            if(bsdf_weight.is_zero()) {
                break;
            }

            r.origin = rec.p;
            r.direction = onb.local_to_world(sample_record.wo);

            throughput = throughput * bsdf_weight;
            eta *= sample_record.eta; // 此处的折射率用于俄罗斯轮盘

            // 俄罗斯轮盘：depth 为 0..4 的弹射一律保留。
            // 从第 6 次弹射起按存活概率 q 随机终止；活下来的路径除以 q，估计保持无偏。
            if(depth >= rr_depth) {
                // 最亮通道代表这条路径的贡献，再乘 eta^2，计入折射界面上的亮度缩放。
                // 上限 0.95 保证至少有 5% 的概率终止，避免全内反射时一直弹到 MAX_DEPTH。
                const auto q = std::min(
                    std::max({throughput.r, throughput.g, throughput.b}) * eta * eta,
                    global_precision{0.95f});

                if(sobol::sampler(sample_index, 4 + depth * 3) >= q) {
                    break;
                }
                throughput /= q; // 存活，则继续传递能量
            }
        }

        return radiance;
    }
}