#pragma once

#include "geometry/triangle_mesh.h"

#include <fstream>
#include <sstream>
#include <algorithm>
#include <string>

namespace LinRT {
    inline int obj_index(int index, int count) {
        return index < 0 ? count + index : index - 1;
    }

    // 包围盒中心移到原点，最长边缩成 1。只改顶点。
    inline bool center_and_unitize(std::vector<vec3>& positions, const std::vector<std::uint32_t>& indices) {
        vec3 bound_min = positions[indices[0]];
        vec3 bound_max = bound_min;
        for(const auto index : indices) {
            const auto& position = positions[index];
            bound_min.x = std::min(bound_min.x, position.x);
            bound_min.y = std::min(bound_min.y, position.y);
            bound_min.z = std::min(bound_min.z, position.z);
            bound_max.x = std::max(bound_max.x, position.x);
            bound_max.y = std::max(bound_max.y, position.y);
            bound_max.z = std::max(bound_max.z, position.z);
        }

        const auto center = (bound_min + bound_max) * 0.5f;
        const auto extent = bound_max - bound_min;
        const auto longest = std::max(extent.x, std::max(extent.y, extent.z));
        if(longest <= 0.0f) {
            return false;
        }

        const auto scale = 1.0f / longest;
        for(auto& position : positions) {
            position = (position - center) * scale;
        }
        return true;
    }

    // 打不开、索引越界、或没有任何三角形时返回空指针。调用处不要换成别的几何体。
    inline std::shared_ptr<triangle_mesh> load_obj(const std::string& path) {
        std::ifstream file(path);
        if(!file) {
            return nullptr;
        }

        std::vector<vec3> positions;
        std::vector<vec3> normals;
        std::vector<std::uint32_t> indices;
        std::vector<vec3> corner_normals;
        std::string line;

        while(std::getline(file, line)) {
            if(!line.empty() && line.back() == '\r') {
                line.pop_back();
            }
            if(line.empty() || line[0] == '#') {
                continue;
            }

            std::istringstream stream(line);
            std::string tag;
            stream >> tag;
            if(tag == "v") {
                float x = 0.0f, y = 0.0f, z = 0.0f;
                stream >> x >> y >> z;
                positions.emplace_back(x, y, z);
            } else if(tag == "vn") {
                float x = 0.0f, y = 0.0f, z = 0.0f;
                stream >> x >> y >> z;
                normals.emplace_back(x, y, z);
            } else if(tag == "f") {
                std::vector<std::uint32_t> face_positions;
                std::vector<int> face_normals;
                std::string corner;
                bool face_ok = true;
                while(stream >> corner) {
                    int position_index = 0;
                    int normal_index = 0;
                    const auto first = corner.find('/');
                    if(first == std::string::npos) {
                        position_index = std::stoi(corner);
                    } else {
                        position_index = std::stoi(corner.substr(0, first));
                        const auto second = corner.find('/', first + 1);
                        if(second != std::string::npos && second + 1 < corner.size()) {
                            normal_index = std::stoi(corner.substr(second + 1));
                        }
                    }

                    const auto position = obj_index(position_index, static_cast<int>(positions.size()));
                    if(position < 0 || position >= static_cast<int>(positions.size())) {
                        face_ok = false;
                        break;
                    }
                    face_positions.push_back(static_cast<std::uint32_t>(position));

                    if(normal_index == 0) {
                        face_normals.push_back(-1);
                    } else {
                        const auto normal = obj_index(normal_index, static_cast<int>(normals.size()));
                        if(normal < 0 || normal >= static_cast<int>(normals.size())) {
                            face_ok = false;
                            break;
                        }
                        face_normals.push_back(normal);
                    }
                }
                if(!face_ok || face_positions.size() < 3) {
                    return nullptr;
                }

                for(size_t i = 1; i + 1 < face_positions.size(); ++i) {
                    const std::uint32_t tri[3] = {face_positions[0], face_positions[i], face_positions[i + 1]};
                    const int tri_normals[3] = {face_normals[0], face_normals[i], face_normals[i + 1]};
                    const auto geometric = cross(positions[tri[1]] - positions[tri[0]], positions[tri[2]] - positions[tri[0]]);
                    if(geometric.length_squared() <= 0.0f) {
                        continue;
                    }
                    const auto unit = geometric / geometric.length();
                    for(int k = 0; k < 3; ++k) {
                        indices.push_back(tri[k]);
                        corner_normals.push_back(tri_normals[k] < 0 ? unit : normals[static_cast<size_t>(tri_normals[k])]);
                    }
                }
            }
        }

        if(indices.empty()) {
            return nullptr;
        }

        if(!center_and_unitize(positions, indices)) {
            return nullptr;
        }

        return std::make_shared<triangle_mesh>(std::move(positions), std::move(indices), std::move(corner_normals));
    }
}