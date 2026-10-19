#include "io/image_write.h"
#include "io/obj_loader.h"
#include "core/scene.h"
#include "core/camera.h"
#include "integrator/path/path_integrator.h"
#include "integrator/photonmapper/photon_integrator.h"
#include "core/hittable_list.h"
#include "geometry/sphere.h"
#include "geometry/rectangle.h"
#include "material/diffuse_bsdf.h"
#include "material/rough_conductor_bsdf.h"
#include "material/dielectric_bsdf.h"
#include "material/rough_dielectric_bsdf.h"
#include "material/plastic_bsdf.h"
#include "material/rough_plastic_bsdf.h"
#include "material/coating_bsdf.h"
#include "accelerator/thread_pool.h"
#include "emitter/area_emitter.h"
#include "emitter/env_emitter.h"
#include "core/sampler.h"

#include <chrono>
#include <iostream>

constexpr int WIDTH = 512;
constexpr int HEIGHT = 512;
constexpr int SAMPLES_PER_PIXEL = 128;

constexpr int TILE = 32;

#if 0
#define USE_PATH_INTEGRATOR
#else
#define USE_PHOTON_MAP_INTEGRATOR
#endif

const auto diffuse = std::make_shared<LinRT::diffuse_bsdf>(LinRT::color(.3f, .6f, .6f));
const auto metal = std::make_shared<LinRT::rough_conductor_bsdf>(LinRT::conductor_ior::copper(), 0.3f);
const auto glass = std::make_shared<LinRT::dielectric_bsdf>(1.5f, 1.0f);
const auto rough_glass = std::make_shared<LinRT::rough_dielectric_bsdf>(1.5f, 1.0f, 0.05f);
const auto plastic = std::make_shared<LinRT::plastic_bsdf>(1.49f, 1.000277f, LinRT::color(1.0f), LinRT::color(0.8f, 0.8f, 0.5f));
const auto rough_plastic = std::make_shared<LinRT::rough_plastic_bsdf>(1.49f, 1.000277f, 0.1f, LinRT::color(1.0f), LinRT::color(0.8f, 0.1f, 0.05f));
const auto coated_copper = std::make_shared<LinRT::coating_bsdf>(
    1.7f, 1.000277f, plastic, LinRT::color(0.3,0.6,0.6));

const auto env_emitter = std::make_shared<LinRT::env_emitter>(
    LinRT::vec3(0.6f, 0.4f, -0.7f), // 太阳方向，内部会归一化
    3.0f,                            // 浑浊度
    0.3f,                            // 地面反照率
    0.05f);                         // 曝光，画面全白就再调小

const auto sphere_emitter = std::make_shared<LinRT::area_emitter>(
    std::make_shared<LinRT::sphere>(LinRT::vec3(1.2f, 0.4f, 0.0f), 0.3f), LinRT::color(4.0f));

int main() {
    LinRT::camera camera;
    camera.origin = LinRT::vec3(.0f,.5f,3.0f);
    camera.update(
        LinRT::vec3(.0f,1.0f,.0f), // up vector
        60.0f,                    // fov
        static_cast<LinRT::global_precision>(WIDTH) 
        / static_cast<LinRT::global_precision>(HEIGHT)); // aspect ratio

    LinRT::hittable_list world;

    auto mesh = LinRT::load_obj("assets/bunny.obj");
    if(!mesh) {
        std::cerr << "obj failed\n";
        return 1;
    }
    mesh->translate(LinRT::vec3(0.0f, 0.0f, 0.0f));
    mesh->set_material(coated_copper);
    world.add(mesh);

    auto plane = std::make_shared<LinRT::sphere>(LinRT::vec3(0.0f, -100.5f, 0.0f), 100.0f);
    plane->set_material(std::make_shared<LinRT::diffuse_bsdf>(LinRT::color(.8f, .8f, .0f)));
    world.add(plane);
    
    auto rect = std::make_shared<LinRT::rectangle>(
        LinRT::vec3(-0.6f, 1.5f, -0.6f),
        LinRT::vec3(1.2f, 0.0f, 0.0f),
        LinRT::vec3(0.0f, 0.0f, 1.2f));
    auto rect_emitter = std::make_shared<LinRT::area_emitter>(rect, LinRT::color(12.0f));        
    world.add(rect);

    std::vector<std::shared_ptr<LinRT::emitter>> emitters{ rect_emitter };
    LinRT::scene render_scene(std::move(world), std::move(emitters), camera);

#ifdef USE_PATH_INTEGRATOR
    const LinRT::path_integrator integrator;
#elif defined(USE_PHOTON_MAP_INTEGRATOR)
    const auto photon_map_build_start = std::chrono::steady_clock::now();
    LinRT::photon_integrator integrator;
    integrator.preprocess(render_scene);
    const auto photon_map_build_end = std::chrono::steady_clock::now();
    const auto photon_map_build_time = std::chrono::duration<double>(photon_map_build_end - photon_map_build_start).count();
    std::cout << "Photon map build time: " << photon_map_build_time << " seconds\n";
#endif

    const auto start = std::chrono::steady_clock::now();

    std::vector<LinRT::color> pixels(WIDTH * HEIGHT);
    const auto total = static_cast<long long>(WIDTH) * HEIGHT;
    std::atomic<long long> done(0);
    std::mutex progress_mutex;
    {
        LinRT::thread_pool pool;
        for(int y0 = 0; y0 < HEIGHT; y0 += TILE) {
            for(int x0 = 0; x0 < WIDTH; x0 += TILE) {
                const int x1 = std::min(x0 + TILE, WIDTH);
                const int y1 = std::min(y0 + TILE, HEIGHT);

                // 块的范围按值拷进任务。场景和相机只读，像素区域互不重叠。
                pool.submit([&, x0, y0, x1, y1]{
                    for(int y = y0; y < y1; ++y) {
                        for(int x = x0; x < x1; ++x) {
                            LinRT::color pixel_color(0.0f);
                            for(int s = 0; s < SAMPLES_PER_PIXEL; ++s) {
                                const auto pixel_index = static_cast<uint32_t>(y * WIDTH + x);
                                LinRT::sampler sampler;
                                sampler.start(pixel_index, static_cast<uint32_t>(s));
                                const auto offset = sampler.next_2d();
                                
                                const auto u = (static_cast<LinRT::global_precision>(x) + offset.x)
                                    / static_cast<LinRT::global_precision>(WIDTH - 1);
                                const auto v = (static_cast<LinRT::global_precision>(y) + offset.y)
                                    / static_cast<LinRT::global_precision>(HEIGHT - 1);
                                
                                auto ray = render_scene.sensor().get_ray(u, v);
                                pixel_color += integrator.radiance(ray, render_scene, sampler);
                            }
                            pixels[static_cast<size_t>(y) * WIDTH + x] =
                                pixel_color / static_cast<LinRT::global_precision>(SAMPLES_PER_PIXEL);
                        }
                    }
                    done.fetch_add(static_cast<long long>(x1 - x0) * (y1 - y0));
                    const auto line = "\r" + std::to_string(done.load()) + "/" + std::to_string(total);
                    std::lock_guard<std::mutex> lock(progress_mutex);
                    std::cout << line << std::flush;
                });
            }
        }
        pool.wait();
    }
    const auto seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    std::cout << "\n" << seconds << " seconds\n";

    LinRT::ImageWriter::save("output.png", WIDTH, HEIGHT, pixels, LinRT::ImageFormat::PNG, true);

    return 0;
}