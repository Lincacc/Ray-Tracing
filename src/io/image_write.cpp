#include "image_write.h"
#include <stb/stb_image_write.h>

#include <iostream>
#include <cmath>
#include <fstream>
#include <algorithm>

namespace LinRT {
    unsigned char to_byte(global_precision value) {
        if(std::isnan(value)) {
            value = 0.0f; // 将 NaN 转换为黑色
        }

        value = std::pow(value, 1.0f / 2.2f); // 伽马校正

        value = std::clamp(value, 0.0f, 0.999f); // 限制在 [0, 1] 范围内

        return static_cast<unsigned char>(value * 256.0f); // 转换为 [0, 255] 范围
    }

    bool ImageWriter::save(const std::string& filepath, 
                            int width, 
                            int height, 
                            const std::vector<color>& pixels, 
                            ImageFormat format, 
                            bool flip_vertically,
                            int quality) {
        if(pixels.size() != static_cast<size_t>(width * height)) {
            std::cerr << "Error: Pixel data size does not match width and height." << std::endl;
            return false;
        }

        // 处理简单纯文本格式图像 PPM，不依赖 STB 图像库
        if(format == ImageFormat::PPM) {
            std::ofstream out(filepath);
            if(!out.is_open()) {
                return false;
            }

            out << "P3\n" << width << " " << height << "\n255\n";
            for(const auto& pixel : pixels) {
                out << static_cast<int>(to_byte(pixel.r)) << " "
                    << static_cast<int>(to_byte(pixel.g)) << " "
                    << static_cast<int>(to_byte(pixel.b)) << "\n";
            }
            out.close();

            return true;
        }

        // 将浮点 RGB 转换成 STB 所需的连续 uint8 字节流 (RGB 3 通道)
        std::vector<unsigned char> byte_data(width * height * 3);
        for(int i = 0; i < width * height; ++i) {
            byte_data[i * 3 + 0] = to_byte(pixels[i].r);
            byte_data[i * 3 + 1] = to_byte(pixels[i].g);
            byte_data[i * 3 + 2] = to_byte(pixels[i].b);
        }

        if(flip_vertically) {
            // 翻转图像数据
            for(int y = 0; y < height / 2; ++y) {
                for(int x = 0; x < width; ++x) {
                    std::swap(byte_data[(y * width + x) * 3 + 0], byte_data[((height - 1 - y) * width + x) * 3 + 0]);
                    std::swap(byte_data[(y * width + x) * 3 + 1], byte_data[((height - 1 - y) * width + x) * 3 + 1]);
                    std::swap(byte_data[(y * width + x) * 3 + 2], byte_data[((height - 1 - y) * width + x) * 3 + 2]);
                }
            }
        }

        int result = 0;
        switch(format) {
            case ImageFormat::PNG:
                result = stbi_write_png(filepath.c_str(), width, height, 3, byte_data.data(), width * 3);
                break;
            case ImageFormat::JPEG:
                result = stbi_write_jpg(filepath.c_str(), width, height, 3, byte_data.data(), quality);
                break;
            case ImageFormat::HDR:
                result = stbi_write_hdr(filepath.c_str(), width, height, 3, reinterpret_cast<float*>(const_cast<color*>(pixels.data())));
                break;
            default:
                std::cerr << "Error: Unsupported image format." << std::endl;
                return false;
        }

        if(!result) {
            std::cerr << "Error: Failed to write image to file." << std::endl;
            return false;
        }
        
        return true;
    }
}