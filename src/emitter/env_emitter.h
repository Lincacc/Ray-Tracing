#pragma once

#include "emitter/emitter.h"
#include "ArHosekSkyModel.h"

namespace LinRT {
    class env_emitter : public emitter {
    public:
        // turbidity：1 到 10，2 很蓝，3 到 5 是常见晴天
        // albedo：地面反照率，0 到 1
        // scale：曝光，用来把模型的线性亮度压进画面
        env_emitter(const vec3& sun_direction, global_precision turbidity, global_precision albedo, global_precision scale)
        : sun_direction(normalize(sun_direction)), scale(scale) {
            const auto elevation = std::asin(std::clamp(this->sun_direction.y, -1.0f, 1.0f));
            state.reset(arhosek_rgb_skymodelstate_alloc_init(turbidity, albedo, elevation));
        }

        bool sample(const vec3& shading_point, const vec2& sample, direct_sample& out) const override {
            const auto z = 1.0f - 2.0f * sample.x;
            const auto r = std::sqrt(std::max(0.0f, 1.0f - z * z));
            const auto phi = 2.0f * M_PI<global_precision> * sample.y;

            out.direction = vec3(
                r * std::cos(phi),
                r * std::sin(phi),
                z);
            out.pdf = 1.0f / (4.0f * M_PI<global_precision>);
            out.radiance = eval(out.direction);

            return true;
        }

        global_precision pdf(const vec3& shading_point, const vec3& direction) const override {
            return 1.0f / (4.0f * M_PI<global_precision>);
        }

        color eval(const vec3& direction) const override {
            const auto dir = normalize(direction);

            // 天顶是 +Y。地平线以下模型没有定义。
            if(dir.y <= 0.0f) {
                return color(0.0f);
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
        global_precision scale = 1.0f;
        std::unique_ptr<ArHosekSkyModelState, state_deleter> state;
    };
}