#pragma once

#include "core/camera.h"
#include "core/hittable_list.h"
#include "emitter/emitter.h"

#include <algorithm>
#include <vector>

namespace LinRT {
    // Mitsuba RGB Spectrum::getLuminance。用来把功率谱收成 DiscretePDF 的标量权重。
    inline global_precision luminance(const color& value) {
        return 0.212671f * value.r + 0.715160f * value.g + 0.072169f * value.b;
    }

    // 场景只保存形状、光源和相机，并提供两种光源抽样。
    // 下一事件对全部光源均匀抽取，和原来的 path 相同。
    // 发射射线只从功率为正的面积光里抽，离散密度正比于 power 的亮度。
    class scene {
    public:
        scene(hittable_list shapes, std::vector<std::shared_ptr<emitter>> emitters, camera sensor)
            : shapes(std::move(shapes)), emitters(std::move(emitters)), sensor_(std::move(sensor)) {
            global_precision sum = 0.0f;
            for(size_t i = 0; i < this->emitters.size(); ++i) {
                const auto weight = luminance(this->emitters[i]->power());
                if(weight <= 0.0f) {
                    continue;
                }
                sum += weight;
                emission_index.push_back(i);
                emission_cdf.push_back(sum);
            }
            emission_power = sum;
        }

        bool hit(const ray& r, global_precision t_min, global_precision t_max, hit_record& rec) const {
            return shapes.hit(r, t_min, t_max, rec);
        }

        const camera& sensor() const { return sensor_; }

        // 飞出场景时的环境辐射。只累加无限远光源的 eval。
        color eval_environment(const vec3& direction) const {
            color emitted(0.0f);
            for(const auto& light : emitters) {
                if(light->is_infinite()) {
                    emitted += light->eval(direction);
                }
            }
            return emitted;
        }

        // 用二维样本的 x 分量先抽中哪一盏灯，余数再交给这盏灯。各灯概率都是 1/n。
        bool sample_emitter_direction(const vec3& shading_point, const vec2& sample, direct_sample& out) const {
            const auto count = static_cast<int>(emitters.size());
            if(count <= 0) {
                return false;
            }

            const auto scaled = sample.x * static_cast<global_precision>(count);
            auto index = static_cast<int>(scaled);
            if(index >= count) {
                index = count - 1;
            }

            const auto resampled = vec2(scaled - static_cast<global_precision>(index), sample.y);
            if(!emitters[static_cast<size_t>(index)]->sample(shading_point, resampled, out)) {
                return false;
            }

            out.emitter_ptr = emitters[static_cast<size_t>(index)].get();
            out.pdf *= 1.0f / static_cast<global_precision>(count);
            return out.pdf > 0.0f;
        }

        // 这个方向被光源策略抽到的总密度。某盏灯抽不到该方向时，它的 pdf 返回 0。
        global_precision pdf_emitter_direction(const vec3& shading_point, const vec3& direction) const {
            const auto count = static_cast<int>(emitters.size());
            if(count <= 0) {
                return 0.0f;
            }

            const auto selection = 1.0f / static_cast<global_precision>(count);
            global_precision sum = 0.0f;
            for(const auto& light : emitters) {
                sum += selection * light->pdf(shading_point, direction);
            }
            return sum;
        }

        // Mitsuba Scene::sampleEmitterRay：先按离散密度选灯，再调用该灯的 sample_ray，权重除以离散密度。
        // 这里的离散密度是 luminance(power) / 总亮度，不是下一事件用的 1/n。
        // 返回后：weight * selection_pdf * pdf_position * pdf_direction = Le * cosθ。
        bool sample_emitter_ray(const vec2& spatial, const vec2& directional, emitter_ray_sample& out) const {
            if(emission_power <= 0.0f || emission_cdf.empty()) {
                return false;
            }

            const auto target = spatial.x * emission_power;
            const auto it = std::upper_bound(emission_cdf.begin(), emission_cdf.end(), target);
            auto bin = static_cast<size_t>(it - emission_cdf.begin());
            if(bin >= emission_cdf.size()) {
                bin = emission_cdf.size() - 1;
            }

            const auto prev = bin == 0 ? 0.0f : emission_cdf[bin - 1];
            const auto span = emission_cdf[bin] - prev;
            if(span <= 0.0f) {
                return false;
            }

            const auto selection_pdf = span / emission_power;
            const auto u = (target - prev) / span;
            const auto& light = emitters[emission_index[bin]];
            if(!light->sample_ray(vec2(u, spatial.y), directional, out)) {
                return false;
            }

            out.weight /= selection_pdf;
            out.selection_pdf = selection_pdf;
            out.emitter_ptr = light.get();
            return true;
        }

    private:
        hittable_list shapes;
        std::vector<std::shared_ptr<emitter>> emitters;
        camera sensor_;
        std::vector<size_t> emission_index;
        std::vector<global_precision> emission_cdf;
        global_precision emission_power = 0.0f;
    };
}
