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
            w = normalize(origin - look_at);
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

        // 把场景中的点投回 get_ray 使用的视口。s、t 与 get_ray 的参数相同，范围 [0, 1]。
        // cos_theta 是「从相机指向该点」与相机朝向 -w 的夹角余弦。
        // 只有针孔（lens_radius == 0）有这个逆变换。光圈不为 0 时射线原点不在一点上。
        bool project(const vec3& point, global_precision& s, global_precision& t, global_precision& cos_theta) const {
            const auto to_point = point - origin;
            const auto dist = to_point.length();
            if(dist <= 0.0f) {
                return false;
            }
            return film_coords(to_point / dist, s, t, cos_theta);
        }

        // Mitsuba PerspectiveCamera::importance。
        // 均匀采样单位距离处的视口矩形，面积是 A'，立体角密度是 1 / (A' * cos³θ)。
        // 方向落在视口外或相机背后时返回 0。
        global_precision pdf_direction(const vec3& world_direction) const {
            const auto length = world_direction.length();
            if(length <= 0.0f) {
                return 0.0f;
            }

            global_precision s = 0.0f;
            global_precision t = 0.0f;
            global_precision cos_theta = 0.0f;
            if(!film_coords(world_direction / length, s, t, cos_theta)) {
                return 0.0f;
            }

            const auto area = unit_rect_area();
            if(area <= 0.0f) {
                return 0.0f;
            }
            return 1.0f / (area * cos_theta * cos_theta * cos_theta);
        }

        // 相机朝向前方、距离原点 1 的那块视口矩形的面积。focus_dist 会在这里约掉。
        global_precision unit_rect_area() const {
            const auto center = lower_left_corner + horizontal * 0.5f + vertical * 0.5f;
            const auto dist2 = (center - origin).length_squared();
            if(dist2 <= 0.0f) {
                return 0.0f;
            }
            return horizontal.length() * vertical.length() / dist2;
        }

    private:
        bool film_coords(const vec3& unit_dir, global_precision& s, global_precision& t,
                        global_precision& cos_theta) const {
            if(lens_radius > 0.0f) {
                return false;
            }

            cos_theta = dot(unit_dir, -w);
            if(cos_theta <= 0.0f) {
                return false;
            }

            const auto center = lower_left_corner + horizontal * 0.5f + vertical * 0.5f;
            const auto focus = (center - origin).length();
            if(focus <= 0.0f) {
                return false;
            }

            // 视口平面过 center，法线是 w。单位方向走到该平面的距离是 focus / cosθ。
            const auto hit = origin + (focus / cos_theta) * unit_dir;
            const auto rel = hit - lower_left_corner;
            const auto hu = dot(horizontal, horizontal);
            const auto hv = dot(vertical, vertical);
            if(hu <= 0.0f || hv <= 0.0f) {
                return false;
            }

            s = dot(rel, horizontal) / hu;
            t = dot(rel, vertical) / hv;
            return s >= 0.0f && s <= 1.0f && t >= 0.0f && t <= 1.0f;
        }

    public:
        vec3 origin;
        vec3 look_at;
        vec3 lower_left_corner;
        vec3 horizontal;
        vec3 vertical;
        vec3 u, v, w;
        global_precision lens_radius = 0.0f;
    };
}