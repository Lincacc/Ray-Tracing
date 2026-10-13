#pragma once

#include "ray.h"
#include "sobol.h"

namespace LinRT {
    class camera {
    public:
        camera() = default;

        /**
         * @brief 更新相机内部参数
         * 
         * @param up            相机的上方向向量
         * @param fov           相机的垂直视场角 (单位: 度)
         * @param aspect_ratio  相机的宽高比 (width / height)
         * @param aperture      相机的光圈大小 (单位: 米)
         * @param focus_dist    相机的焦距 (单位: 米)
         */
        void update(
                    const vec3& up, 
                    global_precision fov, 
                    global_precision aspect_ratio, 
                    global_precision aperture = 0.0f, 
                    global_precision focus_dist = 1.0f) {
            auto h = std::tan(radians(fov) / 2.0f);
            auto viewport_height = 2.0f * h;
            auto viewport_width = aspect_ratio * viewport_height;

            // 构建相机正交基坐标系
            w = normalize(origin - lower_left_corner);
            u = normalize(cross(up, w));
            v = cross(w, u);

            // 计算相机视口参数
            horizontal = focus_dist * viewport_width * u;
            vertical = focus_dist * viewport_height * v;

            // 计算相机视口左下角位置
            lower_left_corner = origin - horizontal / 2.0f - vertical / 2.0f - focus_dist * w;

            // 计算相机光圈半径
            lens_radius = aperture / 2.0f;
        }

        /**
         * @brief 生成一条从相机发出的光线
         * 
         * @param s 水平坐标 (0.0 - 1.0)
         * @param t 垂直坐标 (0.0 - 1.0)
         * @return vec3 光线的方向向量
         */
        ray get_ray(global_precision s, global_precision t) const {
            // 如果光圈半径为 0，则直接生成一条从相机原点发出的光线
            if(lens_radius <= 0.0f) {
                return ray(
                    origin,
                    normalize(lower_left_corner + s * horizontal + t * vertical - origin)
                );
            }

            vec3 random_disk = vec3(sobol::sampler(0, 0), sobol::sampler(0, 1), 0.0f);
            vec3 rd = lens_radius * random_disk;
            vec3 offset = u * rd.x + v * rd.y;

            vec3 ray_origin = origin + offset;
            vec3 ray_target = lower_left_corner + s * horizontal + t * vertical;

            return ray(
                ray_origin,
                normalize(ray_target - ray_origin)
            );
        }
        
    public:
        vec3 origin;
        vec3 lower_left_corner;
        vec3 horizontal;
        vec3 vertical;
        vec3 u, v, w;
        global_precision lens_radius = 0.0f;
    };
}