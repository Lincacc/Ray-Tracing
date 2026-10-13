#include "io/image_write.h"
#include "core/tracer.h"
#include "core/camera.h"
#include "core/hittable_list.h"
#include "geometry/sphere.h"
#include "material/diffuse_bsdf.h"

#include <chrono>
#include <iostream>

constexpr int WIDTH = 512;
constexpr int HEIGHT = 512;
constexpr int SAMPLES_PER_PIXEL = 128;

int main() {
    LinRT::camera camera;
    camera.origin = LinRT::vec3(.0f,.0f,3.0f);
    camera.update(
        LinRT::vec3(.0f,1.0f,.0f), // up vector
        60.0f,                    // fov
        static_cast<LinRT::global_precision>(WIDTH) 
        / static_cast<LinRT::global_precision>(HEIGHT)); // aspect ratio

    LinRT::hittable_list world;
    auto sphere0 = std::make_shared<LinRT::sphere>(LinRT::vec3(0.0f), 0.5f);
    auto diffuse_bsdf = std::make_shared<LinRT::diffuse_bsdf>(LinRT::color(.3f, .6f, .6f));
    sphere0->set_material(diffuse_bsdf);
    world.add(sphere0);

    auto sphere1 = std::make_shared<LinRT::sphere>(LinRT::vec3(0.0f, -100.5f, 0.0f), 100.0f);
    auto diffuse_bsdf1 = std::make_shared<LinRT::diffuse_bsdf>(LinRT::color(.8f, .8f, .0f));
    sphere1->set_material(diffuse_bsdf1);
    world.add(sphere1);
    
    auto env_emitter = std::make_shared<LinRT::env_emitter>(
        LinRT::vec3(0.6f, 0.4f, -0.7f), // 太阳方向，内部会归一化
        3.0f,                            // 浑浊度
        0.3f,                            // 地面反照率
        0.05f);                         // 曝光，画面全白就再调小

    const auto start = std::chrono::steady_clock::now();

    std::vector<LinRT::color> pixels(WIDTH * HEIGHT);
    for(int y = 0; y < HEIGHT; ++y) {
        for(int x = 0; x < WIDTH; ++x) {

            LinRT::color pixel_color(0.0f);
            LinRT::hit_record rec;
            for(int s = 0; s < SAMPLES_PER_PIXEL; ++s) {
                auto u_offset = LinRT::sobol::sampler(s, 0);
                auto v_offset = LinRT::sobol::sampler(s, 1);

                auto u = (static_cast<LinRT::global_precision>(x) + u_offset) / static_cast<LinRT::global_precision>(WIDTH - 1);
                auto v = (static_cast<LinRT::global_precision>(y) + v_offset) / static_cast<LinRT::global_precision>(HEIGHT - 1);

                auto ray = camera.get_ray(u, v);
                pixel_color += LinRT::trace(ray, world, env_emitter, static_cast<uint32_t>(s));
            }
            pixel_color /= static_cast<LinRT::global_precision>(SAMPLES_PER_PIXEL);

            pixels[y * WIDTH + x] = pixel_color;

            const auto done = static_cast<long long>(y + 1) * WIDTH;
            const auto total = static_cast<long long>(WIDTH) * HEIGHT;
            const auto line = "\r" + std::to_string(done) + "/" + std::to_string(total);
            std::cout << line << std::flush;
        }
    }

    const auto seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    std::cout << "\n" << seconds << " seconds\n";

    LinRT::ImageWriter::save("output.png", WIDTH, HEIGHT, pixels, LinRT::ImageFormat::PNG, true);

    return 0;
}