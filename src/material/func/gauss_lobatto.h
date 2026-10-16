#pragma once

#include "core/core.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>

namespace LinRT {
    // Mitsuba src/libcore/quad.cpp 的 GaussLobattoIntegrator。
    // SmoothPlastic::configure 用它算 fresnelDiffuseReflectance(eta, false)：
    // 最多 1024 次求值，绝对误差 0，相对误差 1e-5，并打开收敛估计。
    // 停机条件 dist == acc 来自 QuantLib：差值已经落进 acc 的一个 ulp。
    class gauss_lobatto_integrator {
    public:
        gauss_lobatto_integrator(std::size_t max_evals, global_precision abs_error, global_precision rel_error)
            : max_evals(max_evals)
            , abs_error(abs_error)
            , rel_error(rel_error) {}

        template<typename Integrand>
        global_precision integrate(const Integrand& f, global_precision a, global_precision b) const {
            auto factor = global_precision{1};
            std::size_t evals = 0;
            if(a == b) {
                return 0.0f;
            }
            if(b < a) {
                std::swap(a, b);
                factor = -1.0f;
            }

            const auto tolerance = absolute_tolerance(f, a, b, evals);
            evals += 2;
            return factor * adaptive_step(f, a, b, f(a), f(b), tolerance, evals);
        }

    private:
        template<typename Integrand>
        global_precision absolute_tolerance(const Integrand& f, global_precision a, global_precision b,
                                            std::size_t& evals) const {
            const auto m = (a + b) * 0.5f;
            const auto h = (b - a) * 0.5f;
            const auto y1 = f(a);
            const auto y3 = f(m - alpha * h);
            const auto y5 = f(m - beta * h);
            const auto y7 = f(m);
            const auto y9 = f(m + beta * h);
            const auto y11 = f(m + alpha * h);
            const auto y13 = f(b);

            const auto acc = h * (static_cast<global_precision>(0.0158271919734801831) * (y1 + y13)
                + static_cast<global_precision>(0.0942738402188500455) * (f(m - x1 * h) + f(m + x1 * h))
                + static_cast<global_precision>(0.1550719873365853963) * (y3 + y11)
                + static_cast<global_precision>(0.1888215739601824544) * (f(m - x2 * h) + f(m + x2 * h))
                + static_cast<global_precision>(0.1997734052268585268) * (y5 + y9)
                + static_cast<global_precision>(0.2249264653333395270) * (f(m - x3 * h) + f(m + x3 * h))
                + static_cast<global_precision>(0.2426110719014077338) * y7);
            evals += 13;

            auto r = global_precision{1};
            const auto integral2 = (h / 6.0f) * (y1 + y13 + 5.0f * (y5 + y9));
            const auto integral1 = (h / 1470.0f) * (77.0f * (y1 + y13) + 432.0f * (y3 + y11)
                + 625.0f * (y5 + y9) + 672.0f * y7);
            if(std::abs(integral2 - acc) != 0.0f) {
                r = std::abs(integral1 - acc) / std::abs(integral2 - acc);
            }
            if(r == 0.0f || r > 1.0f) {
                r = 1.0f;
            }

            const auto epsilon = std::numeric_limits<global_precision>::epsilon();
            auto result = std::numeric_limits<global_precision>::infinity();
            if(rel_error != 0.0f && acc != 0.0f) {
                result = acc * std::max(rel_error, epsilon) / (r * epsilon);
            }
            if(abs_error != 0.0f) {
                result = std::min(result, abs_error / (r * epsilon));
            }
            return result;
        }

        template<typename Integrand>
        global_precision adaptive_step(const Integrand& f, global_precision a, global_precision b,
                                       global_precision fa, global_precision fb, global_precision acc,
                                       std::size_t& evals) const {
            const auto h = (b - a) * 0.5f;
            const auto m = (a + b) * 0.5f;
            const auto mll = m - alpha * h;
            const auto ml = m - beta * h;
            const auto mr = m + beta * h;
            const auto mrr = m + alpha * h;

            const auto fmll = f(mll);
            const auto fml = f(ml);
            const auto fm = f(m);
            const auto fmr = f(mr);
            const auto fmrr = f(mrr);

            // 4 点 Lobatto，以及补上中点和两个 alpha 节点之后的 7 点 Kronrod。
            const auto integral2 = (h / 6.0f) * (fa + fb + 5.0f * (fml + fmr));
            const auto integral1 = (h / 1470.0f) * (77.0f * (fa + fb) + 432.0f * (fmll + fmrr)
                + 625.0f * (fml + fmr) + 672.0f * fm);
            evals += 5;
            if(evals >= max_evals) {
                return integral1;
            }

            const auto dist = acc + (integral1 - integral2);
            if(dist == acc || mll <= a || b <= mrr) {
                return integral1;
            }

            return adaptive_step(f, a, mll, fa, fmll, acc, evals)
                + adaptive_step(f, mll, ml, fmll, fml, acc, evals)
                + adaptive_step(f, ml, m, fml, fm, acc, evals)
                + adaptive_step(f, m, mr, fm, fmr, acc, evals)
                + adaptive_step(f, mr, mrr, fmr, fmrr, acc, evals)
                + adaptive_step(f, mrr, b, fmrr, fb, acc, evals);
        }

        std::size_t max_evals;
        global_precision abs_error;
        global_precision rel_error;

        // 与 Mitsuba 一样，先写成 double，再收成 Float。直接写 float 字面量会在另一处舍入。
        static inline const global_precision alpha = static_cast<global_precision>(std::sqrt(2.0 / 3.0));
        static inline const global_precision beta = static_cast<global_precision>(1.0 / std::sqrt(5.0));
        static inline const global_precision x1 = static_cast<global_precision>(0.94288241569547971906);
        static inline const global_precision x2 = static_cast<global_precision>(0.64185334234578130578);
        static inline const global_precision x3 = static_cast<global_precision>(0.23638319966214988028);
    };
}
