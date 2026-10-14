#pragma once

#include "geometry/hittable.h"

namespace LinRT {
    class hittable_list{
    public:
        hittable_list() = default;
        hittable_list(std::shared_ptr<hittable> object) { add(object); }

        void add(std::shared_ptr<hittable> object) { objects.push_back(object); }
        void clear() { objects.clear(); }

        bool hit(const ray& r, global_precision t_min, global_precision t_max, hit_record& rec) const {
            hit_record temp_rec;
            bool hit_anything = false;
            // 用于记录当前最近的交点距离
            global_precision closest_so_far = t_max;

            for(const auto& object : objects) {
                if(object->hit(r, t_min, closest_so_far, temp_rec)) {
                    hit_anything = true;
                    closest_so_far = temp_rec.t;
                    rec = temp_rec;
                }
            }

            return hit_anything;
        }

    private:
        std::vector<std::shared_ptr<hittable>> objects;
    };
}