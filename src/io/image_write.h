#pragma once

#include "core/color.h"

#include <string>
#include <vector>

namespace LinRT {
    enum class ImageFormat {
        PNG,
        JPEG,
        PPM,
        HDR
    };

    class ImageWriter {
    public:
        ImageWriter() = default;

        /**
         * @brief 将光追渲染得到的 32 位浮点 RGB 颜色缓冲区导出为图像文件
         * 
         * @param filepath  导出文件路径 (例如 "output.png")
         * @param width     图像宽度
         * @param height    图像高度
         * @param pixels    像素颜色数组（大小为 width * height，按从左到右、从上到下排布）
         * @param format    导出格式 (默认 PNG)
         * @param quality   JPG 压缩质量 (1-100)，仅对 JPG 生效
         * @return true     保存成功
         * @return false    保存失败
         */
        static bool save(const std::string& filepath, 
                        int width, 
                        int height, 
                        const std::vector<color>& pixels, 
                        ImageFormat format = ImageFormat::PNG, 
                        bool flip_vertically = false,
                        int quality = 95);
        };
}