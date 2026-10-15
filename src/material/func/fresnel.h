#pragma once

#include "core/color.h"

#include <algorithm>
#include <complex>

namespace LinRT {
    struct conductor_ior {
        color eta;
        color k;
        
        static conductor_ior copper() {
            return {
                color(0.200438f, 0.924033f, 1.102212f),
                color(3.912948f, 2.452847f, 2.142188f)
            };
        }

        static conductor_ior gold() {
            return {
                color(0.143119f, 0.374957f, 1.442478f),
                color(3.983160f, 2.385721f, 1.603215f)
            };
        }

        static conductor_ior aluminum() {
            return {
                color(1.345600f, 0.965210f, 0.617220f),
                color(7.474600f, 6.399500f, 5.303100f)
            };
        }
        static conductor_ior silver() {
            return {
                color(0.155265f, 0.116723f, 0.138381f),
                color(4.828350f, 3.122245f, 2.146900f)
            };
        }
    };

    inline global_precision fresnel_conductor(global_precision cos_theta_i, global_precision eta, global_precision k) {
        cos_theta_i = std::clamp(cos_theta_i, 0.0f, 1.0f);

        const std::complex<global_precision> ior(eta, k);
        const auto sin2_theta_i = 1.0f - cos_theta_i * cos_theta_i;
        const auto sin2_theta_t = sin2_theta_i / (ior * ior);
        const auto cos_theta_t = std::sqrt(std::complex<global_precision>(1.0f) - sin2_theta_t);

        const auto r_parallel = (ior * cos_theta_i - cos_theta_t) / (ior * cos_theta_i + cos_theta_t);
        const auto r_perpendicular = (cos_theta_i - ior * cos_theta_t) /(cos_theta_i + ior * cos_theta_t);

        // std::norm 是模的平方，振幅反射率变成能量反射率
        return 0.5f * (std::norm(r_parallel) + std::norm(r_perpendicular));
    }

    inline color fresnel_conductor(global_precision cos_theta_i, const color& eta, const color& k) {
        return color(
            fresnel_conductor(cos_theta_i, eta.r, k.r),
            fresnel_conductor(cos_theta_i, eta.g, k.g),
            fresnel_conductor(cos_theta_i, eta.b, k.b)
        );
    }
}