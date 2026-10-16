#pragma once

#include "core/scene.h"
#include "core/onb.h"
#include "core/sampler.h"

#include <algorithm>
#include <cmath>

/*
    从相机出发，沿途做直接光抽样，再做 BSDF 采样。
*/

namespace LinRT {
    // Mitsuba 的 path。弹射、直接光和幂启发式与原来的 trace 相同，选灯改走 scene。
    class path_integrator {
    public:
        color radiance(ray& r, const scene& world, sampler& sampler) const {
            constexpr int rr_depth = 5;
            global_precision eta = 1.0f;
            color throughput(1.0f);
            color result(0.0f);

            // 当前这根射线是不是上一次材质采样产生的，以及那次采样的立体角密度。
            // 相机射线没有这个密度：直接看到天空时，亮度全部由这条射线贡献。
            bool have_solid_angle_pdf = false;
            global_precision solid_angle_pdf = 0.0f;

            for(int depth = 0; depth < max_depth; ++depth) {
                hit_record rec;

                if(!world.hit(r, 0.001f, INFINITY, rec)) {
                    const auto Le = world.eval_environment(r.direction);

                    if(!have_solid_angle_pdf) {
                        result += throughput * Le;
                    } else {
                        // 这条方向如果也落在太阳圆锥里，光源策略同样能抽到它。
                        // 圆锥外 light_pdf 为 0，权重为 1，平滑天光仍只算材质这一路。
                        const auto light_pdf = world.pdf_emitter_direction(rec.p, r.direction);
                        const auto weight = power_heuristic(solid_angle_pdf, light_pdf);
                        result += throughput * weight * Le;
                    }
                    break;
                }

                if(rec.emitter_ptr && rec.front_face) {
                    const auto outgoing = -normalize(r.direction);
                    const auto Le = rec.emitter_ptr->emitted(rec.normal, outgoing);
                    if(!Le.is_zero()) {
                        if(!have_solid_angle_pdf) {
                            result += throughput * Le;
                        } else {
                            const auto light_pdf = world.pdf_emitter_direction(r.origin, r.direction);
                            const auto weight = power_heuristic(solid_angle_pdf, light_pdf);
                            result += throughput * weight * Le;
                        }
                    }
                }

                if(!rec.bsdf_ptr) {
                    break;
                }

                onb basis;
                basis.build_from_normal(rec.normal);

                bsdf_sample_record sample_record;
                sample_record.p = rec.p;
                sample_record.n = normalize(rec.normal);
                // mitsuba: wi 指向上一跳的来源方向
                sample_record.wi = basis.world_to_local(-normalize(r.direction));
                sample_record.front_face = rec.front_face;

                bsdf_context ctx;
                ctx.mode = bsdf_context::transport_mode::RADIANCE;

                // 每个弹射固定 5 维：光源两个，材质两个，俄罗斯轮盘一个。
                // 和 Mitsuba path 一样，没用到也要取走，否则下一次弹射的维度会错位。
                const auto light_uv = sampler.next_2d();
                ctx.lobe_sample = sampler.next_1d();
                const auto bsdf_uv = sampler.next_2d();
                const auto rr_u = sampler.next_1d();

                // 直接光：只对有漫反射或光泽波瓣的材质做。纯镜面没有立体角密度，不能和光源 MIS。
                constexpr unsigned int smooth = DIFFUSE | GLOSSY;
                if((rec.bsdf_ptr->flags() & smooth) != 0) {
                    direct_sample light;

                    if(world.sample_emitter_direction(rec.p, light_uv, light)) {
                        // 起点沿法线挪开之后，必须按新起点重算到灯面上那一点的距离。
                        // 继续用原来的 light.distance，近处的射线会打进灯球自己。
                        const auto shadow_origin = rec.p + rec.normal * 1e-4f;
                        vec3 shadow_direction = light.direction;
                        auto tmax = INFINITY;
                        bool have_target = true;
                        if(!std::isinf(light.distance)) {
                            const auto to_light = light.point - shadow_origin;
                            const auto dist = to_light.length();
                            if(dist < 1e-6f) {
                                have_target = false;
                            } else {
                                shadow_direction = to_light / dist;
                                tmax = dist * (1.0f - 1e-4f);
                            }
                        }

                        if(have_target) {
                            sample_record.wo = basis.world_to_local(shadow_direction);
                            const auto f_cos = rec.bsdf_ptr->eval(ctx, sample_record);
                            if(!f_cos.is_zero()) {
                                hit_record occluder;
                                const ray shadow(shadow_origin, shadow_direction);
                                const bool hit_something = world.hit(shadow, 0.001f, tmax, occluder);
                                // 最近的交点就是目标灯：灯可见。中间先碰到别的物体才算挡住。
                                const bool blocked = hit_something && occluder.emitter_ptr != light.emitter_ptr;
                                if(!blocked) {
                                    const auto bsdf_pdf = rec.bsdf_ptr->pdf(ctx, sample_record);
                                    const auto light_pdf = world.pdf_emitter_direction(rec.p, shadow_direction);
                                    if(light_pdf > 0.0f) {
                                        const auto weight = power_heuristic(light_pdf, bsdf_pdf);
                                        result += throughput * f_cos * light.radiance * (weight / light_pdf);
                                    }
                                }
                            }
                        }
                    }
                }

                // 最后一次弹射只做直接光，不再往下散射。
                if(depth == max_depth - 1) {
                    break;
                }

                auto bsdf_weight = rec.bsdf_ptr->sample(ctx, sample_record, bsdf_uv);
                if(bsdf_weight.is_zero()) {
                    break;
                }

                r.origin = rec.p;
                r.direction = basis.local_to_world(sample_record.wo);

                throughput *= bsdf_weight;
                eta *= sample_record.eta; // 此处的折射率用于俄罗斯轮盘

                // 纯镜面（带 SPECULAR、不带 DIFFUSE/GLOSSY）的密度不是立体角密度，下一次打到天空时不和光源比。
                const bool delta = (sample_record.sampled_type & SPECULAR) != 0
                    && (sample_record.sampled_type & (DIFFUSE | GLOSSY)) == 0;
                have_solid_angle_pdf = !delta && sample_record.pdf > 0.0f;
                solid_angle_pdf = sample_record.pdf;

                // 俄罗斯轮盘：depth 为 0..4 的弹射一律保留。
                // 从第 6 次弹射起按存活概率 q 随机终止；活下来的路径除以 q，估计保持无偏。
                if(depth >= rr_depth) {
                    // 最亮通道代表这条路径的贡献，再乘 eta^2，计入折射界面上的亮度缩放。
                    // 上限 0.95 保证至少有 5% 的概率终止，避免全内反射时一直弹到 max_depth。
                    const auto q = std::min(
                        std::max({throughput.r, throughput.g, throughput.b}) * eta * eta,
                        global_precision{0.95f});

                    if(rr_u >= q) {
                        break;
                    }
                    throughput /= q; // 存活，则继续传递能量
                }
            }

            return result;
        }

    private:
        static constexpr int max_depth = 50;

        // 幂启发式，指数 β = 2，是 Mitsuba path 的默认 MIS。
        // pdf_strategy 是「正在累加的那条策略」的立体角密度，pdf_other 是另一条的。
        static global_precision power_heuristic(global_precision pdf_strategy, global_precision pdf_other) {
            const auto a = pdf_strategy * pdf_strategy;
            const auto b = pdf_other * pdf_other;
            const auto sum = a + b;
            return sum > 0.0f ? a / sum : 0.0f;
        }
    };
}
