#pragma once

#include "emitter/emitter.h"
#include "ArHosekSkyModel.h"
#include "core/onb.h"

namespace LinRT {
    class env_emitter : public emitter {
    public:
        // turbidity：1 到 10，2 很蓝，3 到 5 是常见晴天
        // albedo：地面反照率，0 到 1
        // scale：曝光，用来把模型的线性亮度压进画面
        env_emitter(const vec3& sun_direction, global_precision turbidity, global_precision albedo, global_precision scale)
        : sun_direction(normalize(sun_direction)), scale(scale) {
            local_frame.build_from_normal(this->sun_direction);
            // 太阳在地平线以下时，模型的高亮区不在可见半球里，不再朝它抽样。
            sun_active = this->sun_direction.y > 0.0f;

            // Hosek–Wilkie 使用的太阳视半径是 0.255 度。用 double 计算
            // 1-cos(θ)，避免这个角太小，float 把圆锥立体角吃成 0。
            const double radius = 0.255 * M_PI<double> / 180.0;
            const double cosine = std::cos(radius);
            cos_cutoff = static_cast<global_precision>(cosine);
            // 均匀圆锥的立体角是 2π(1-cos θ)，密度是它的倒数，圆锥内处处相同。
            cone_pdf = static_cast<global_precision>(
                1.0 / (2.0 * M_PI<double> * (1.0 - cosine)));

            const auto elevation = std::asin(std::clamp(this->sun_direction.y, -1.0f, 1.0f));
            state.reset(arhosek_rgb_skymodelstate_alloc_init(turbidity, albedo, elevation));
        }

        // 无限远光源的密度与着色点位置无关。area light 以后才会用这个点。
        bool sample(const vec3& shading_point, const vec2& sample, direct_sample& out) const override {
            (void)shading_point;
            if(!sun_active) {
                return false;
            }

            // Mitsuba warp::square_to_uniform_cone：sample.x = 0 在太阳中心，1 在圆锥边缘。
            const auto cos_theta = (1.0f - sample.x) + sample.x * cos_cutoff;
            const auto sin_theta = std::sqrt(std::max(0.0f, 1.0f - cos_theta * cos_theta));
            const auto phi = 2.0f * M_PI<global_precision> * sample.y;
            const vec3 local(
                std::cos(phi) * sin_theta,
                std::sin(phi) * sin_theta,
                cos_theta);

            out.direction = local_frame.local_to_world(local);
            out.pdf = cone_pdf;
            out.radiance = eval(out.direction);

            return true;
        }

        global_precision pdf(const vec3& shading_point, const vec3& direction) const override {
            (void)shading_point;
            if(!sun_active) {
                return 0.0f;
            }

            const auto dir = normalize(direction);
            // 圆锥边界上，点积可能比 cos_cutoff 小一个舍入误差。留一点余量，
            // 让 sample 刚抽到的方向再问 pdf 时不会突然变成 0。   
            if(dot(dir, this->sun_direction) + 1e-5f < cos_cutoff) {
                return 0.0f;
            }
            
            return cone_pdf;
        }

        color eval(const vec3& direction) const override {
            const auto dir = normalize(direction);

            // 天顶是 +Y。地平线以下模型没有定义。
            if(dir.y <= 0.0f) {
                return color(1.0f);
            }

            // θ：观察方向离天顶的夹角。γ：观察方向和太阳方向的夹角。
            const auto theta = std::acos(std::clamp(dir.y, 0.0f, 1.0f));
            const auto gamma = std::acos(std::clamp(dot(dir, sun_direction), -1.0f, 1.0f));

            return color(
                static_cast<global_precision>(arhosek_tristim_skymodel_radiance(state.get(), theta, gamma, 0)),
                static_cast<global_precision>(arhosek_tristim_skymodel_radiance(state.get(), theta, gamma, 1)),
                static_cast<global_precision>(arhosek_tristim_skymodel_radiance(state.get(), theta, gamma, 2))
            ) * scale;
        }

        bool is_infinite() const override {
            return true;
        }

    private:
        struct state_deleter {
            void operator()(ArHosekSkyModelState* p) const { arhosekskymodelstate_free(p); }
        };

        vec3 sun_direction;
        onb local_frame;
        global_precision scale = 1.0f;
        global_precision cos_cutoff = 1.0f;
        global_precision cone_pdf = 0.0f;
        bool sun_active = false;
        std::unique_ptr<ArHosekSkyModelState, state_deleter> state;
    };
}