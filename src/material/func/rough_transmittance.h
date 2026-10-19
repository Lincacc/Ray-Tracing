#pragma once

#include "material/func/catmull_rom.h"

#include <cstdint>
#include <cstring>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace LinRT {
    // Mitsuba src/bsdfs/rtrans.h。读 ggx.dat，再按 setEta / setAlpha 降到 1D。
    class rough_transmittance {
    public:
        // 路径默认和 Mitsuba 文件解析器相同：data/microfacet/ggx.dat
        rough_transmittance(global_precision eta, global_precision alpha,
                            const std::string& path = "third_party/data/ggx.dat") {
            load(path);
            check_eta(eta);
            check_alpha(alpha);
            set_eta(eta);
            set_alpha(alpha);
        }

        global_precision eval(global_precision cos_theta) const {
            // eta、alpha 都固定之后只查 cosθ。背面为 0。
            if(!(cos_theta >= 0.0f)) {
                return 0.0f;
            }
            const auto warped = std::pow(cos_theta, global_precision{0.25f});
            const auto result = eval_cubic_interp_1d(warped, trans.data(), theta_samples, 0.0f, 1.0f);
            return std::clamp(result, 0.0f, 1.0f);
        }

        global_precision eval_diffuse() const {
            return std::clamp(diff_trans[0], 0.0f, 1.0f);
        }

    private:
        static std::uint64_t read_u64_le(std::ifstream& in) {
            unsigned char b[8];
            in.read(reinterpret_cast<char*>(b), 8);
            return static_cast<std::uint64_t>(b[0])
                | (static_cast<std::uint64_t>(b[1]) << 8)
                | (static_cast<std::uint64_t>(b[2]) << 16)
                | (static_cast<std::uint64_t>(b[3]) << 24)
                | (static_cast<std::uint64_t>(b[4]) << 32)
                | (static_cast<std::uint64_t>(b[5]) << 40)
                | (static_cast<std::uint64_t>(b[6]) << 48)
                | (static_cast<std::uint64_t>(b[7]) << 56);
        }

        // Stream::readSingle：文件里是 little-endian float32。
        static float read_f32_le(std::ifstream& in) {
            unsigned char b[4];
            in.read(reinterpret_cast<char*>(b), 4);
            const auto bits = static_cast<std::uint32_t>(b[0])
                | (static_cast<std::uint32_t>(b[1]) << 8)
                | (static_cast<std::uint32_t>(b[2]) << 16)
                | (static_cast<std::uint32_t>(b[3]) << 24);
            float value;
            std::memcpy(&value, &bits, 4);
            return value;
        }

        void load(const std::string& path) {
            std::ifstream in(path, std::ios::binary);
            if(!in) {
                throw std::runtime_error("rough_transmittance: missing " + path);
            }

            char header[17];
            in.read(header, 17);
            if(!in || std::memcmp(header, "MTS_TRANSMITTANCE", 17) != 0) {
                throw std::runtime_error("rough_transmittance: invalid header in " + path);
            }

            // Stream::readSize = little-endian uint64。
            eta_samples = static_cast<std::size_t>(read_u64_le(in));
            alpha_samples = static_cast<std::size_t>(read_u64_le(in));
            theta_samples = static_cast<std::size_t>(read_u64_le(in));
            eta_min = read_f32_le(in);
            eta_max = read_f32_le(in);
            alpha_min = read_f32_le(in);
            alpha_max = read_f32_le(in);

            // 两块：eta≥1 和 eta<1。每块按 (eta, alpha) 存 theta 个方向透过率，再跟一个漫透过率。
            const auto trans_size = 2 * eta_samples * alpha_samples * theta_samples;
            const auto diff_size = 2 * eta_samples * alpha_samples;
            std::vector<float> packed(trans_size + diff_size);
            for(auto& sample : packed) {
                sample = read_f32_le(in);
            }
            if(!in) {
                throw std::runtime_error("rough_transmittance: truncated " + path);
            }

            trans.resize(trans_size);
            diff_trans.resize(diff_size);
            std::size_t src = 0;
            std::size_t dst_t = 0;
            std::size_t dst_d = 0;
            for(std::size_t i = 0; i < 2 * eta_samples; ++i) {
                for(std::size_t j = 0; j < alpha_samples; ++j) {
                    for(std::size_t k = 0; k < theta_samples; ++k) {
                        trans[dst_t++] = packed[src++];
                    }
                    diff_trans[dst_d++] = packed[src++];
                }
            }
        }

        void check_eta(global_precision eta) const {
            if(eta < 1.0f) {
                eta = 1.0f / eta;
            }
            if(eta < eta_min || eta > eta_max) {
                throw std::runtime_error("rough_transmittance: eta outside table range");
            }
        }

        void check_alpha(global_precision alpha) const {
            if(alpha < alpha_min || alpha > alpha_max) {
                throw std::runtime_error("rough_transmittance: alpha outside table range");
            }
        }

        void set_eta(global_precision eta) {
            auto* trans_src = trans.data();
            auto* diff_src = diff_trans.data();
            if(eta < 1.0f) {
                trans_src += eta_samples * alpha_samples * theta_samples;
                diff_src += eta_samples * alpha_samples;
                eta = 1.0f / eta;
            }
            if(eta < eta_min) {
                eta = eta_min;
            }

            const auto warped_eta = std::pow((eta - eta_min) / (eta_max - eta_min), global_precision{0.25f});
            const auto d_alpha = 1.0f / static_cast<global_precision>(alpha_samples - 1);
            const auto d_theta = 1.0f / static_cast<global_precision>(theta_samples - 1);
            const std::size_t size3[3] = { theta_samples, alpha_samples, eta_samples };
            const std::size_t size2[2] = { alpha_samples, eta_samples };
            const global_precision zero3[3] = { 0.0f, 0.0f, 0.0f };
            const global_precision one3[3] = { 1.0f, 1.0f, 1.0f };
            const global_precision zero2[2] = { 0.0f, 0.0f };
            const global_precision one2[2] = { 1.0f, 1.0f };

            std::vector<global_precision> sliced_t(alpha_samples * theta_samples);
            std::vector<global_precision> sliced_d(alpha_samples);
            for(std::size_t i = 0; i < alpha_samples; ++i) {
                for(std::size_t j = 0; j < theta_samples; ++j) {
                    const global_precision p[3] = {
                        static_cast<global_precision>(j) * d_theta,
                        static_cast<global_precision>(i) * d_alpha,
                        warped_eta
                    };
                    sliced_t[i * theta_samples + j] = eval_cubic_interp_3d(p, trans_src, size3, zero3, one3);
                }
                const global_precision q[2] = { static_cast<global_precision>(i) * d_alpha, warped_eta };
                sliced_d[i] = eval_cubic_interp_2d(q, diff_src, size2, zero2, one2);
            }
            trans = std::move(sliced_t);
            diff_trans = std::move(sliced_d);
        }

        void set_alpha(global_precision alpha) {
            const auto warped_alpha = std::pow((alpha - alpha_min) / (alpha_max - alpha_min),
                                               global_precision{0.25f});
            const auto d_theta = 1.0f / static_cast<global_precision>(theta_samples - 1);
            const std::size_t size2[2] = { theta_samples, alpha_samples };
            const global_precision zero2[2] = { 0.0f, 0.0f };
            const global_precision one2[2] = { 1.0f, 1.0f };

            std::vector<global_precision> sliced_t(theta_samples);
            for(std::size_t i = 0; i < theta_samples; ++i) {
                const global_precision p[2] = { static_cast<global_precision>(i) * d_theta, warped_alpha };
                sliced_t[i] = eval_cubic_interp_2d(p, trans.data(), size2, zero2, one2);
            }
            const auto sliced_d = eval_cubic_interp_1d(warped_alpha, diff_trans.data(), alpha_samples, 0.0f, 1.0f);
            trans = std::move(sliced_t);
            diff_trans.assign(1, sliced_d);
        }

        std::size_t eta_samples = 0;
        std::size_t alpha_samples = 0;
        std::size_t theta_samples = 0;
        global_precision eta_min = 0.0f;
        global_precision eta_max = 0.0f;
        global_precision alpha_min = 0.0f;
        global_precision alpha_max = 0.0f;
        std::vector<global_precision> trans;
        std::vector<global_precision> diff_trans;
    };
}