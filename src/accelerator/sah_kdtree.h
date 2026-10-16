#pragma once

#include "accelerator/aabb.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <utility>
#include <vector>

namespace LinRT {
    // Mitsuba 0.6 加速三角形用的是 SAH kd-tree，不是每个节点带两个子包围盒的二叉树 BVH。
    // 对应 GenericKDTree（gkdtree.h）和 SAHKDTree3D::rayIntersectPBRT（sahkdtree3.h）。
    //
    // 内部节点只保存一个分裂轴和一张分裂平面。平面把父体素切成左右两个长方体。
    // 构建时对每条轴收集图元包围盒的端点事件，从左往右扫描，用表面积启发式选代价最小的平面：
    //   C = 15 + 20 * (P_L * N_L + P_R * N_R)
    // P_L、P_R 是左右子体素表面积除以父体素表面积，和 SurfaceAreaHeuristic3 相同。
    // 某一侧没有图元时，代价再乘 0.9（emptySpaceBonus）。
    // 落在平面上的图元分别试算放左、放右，取得较小的一侧。
    // 叶子代价是 N * 20。划分不比叶子便宜，并且已经明显更差或连续 3 次较差时，就做成叶子。
    // 子树返回后若实际代价仍不低于叶子，撤掉这次分裂（retract）。
    // 三角形数不超过 6，或深度达到 min(48, 8 + 1.3 log2 N)，停止。
    //
    // clip 取 false：跨在平面上的三角形同时放进左右孩子。
    // Mitsuba 默认还会把三角形裁进子体素再生成事件，那一步这次不写。
    // 三角形数超过 65536 时 Mitsuba 会先做 128 箱的 min-max binning，这次一律用上面的事件扫描。
    //
    // 遍历与 rayIntersectPBRT 相同：原点在平面左侧就先进入左孩子；
    // 分裂平面落在当前射线区间内时，才把另一侧压栈。
    class sah_kdtree {
    public:
        void build(const std::vector<aabb>& primitive_bounds) {
            nodes.clear();
            ordered.clear();
            bounds = primitive_bounds;
            const auto count = static_cast<int>(bounds.size());
            if(count <= 0) {
                return;
            }

            aabb root;
            for(const auto& box : bounds) {
                root.extend(box);
            }

            std::vector<std::uint32_t> prims(static_cast<size_t>(count));
            for(int i = 0; i < count; ++i) {
                prims[static_cast<size_t>(i)] = static_cast<std::uint32_t>(i);
            }

            const auto depth_limit = std::min(48, static_cast<int>(8.0f + 1.3f * std::log2(static_cast<float>(count))));
            max_depth = std::max(depth_limit, 1);
            build_node(prims, root, 1, 0);
        }

        // test(三角形下标, 当前体素区间近端, 当前体素区间远端)。
        // 命中且比 closest 更近时，调用者要把 closest 写成该交点距离。
        template <typename Leaf>
        void intersect(const ray& r, global_precision t_min, global_precision& closest, Leaf&& test) const {
            if(nodes.empty()) {
                return;
            }

            struct item {
                int node = 0;
                global_precision mint = 0.0f;
                global_precision maxt = 0.0f;
            };

            std::vector<item> stack;
            auto mint = t_min;
            auto maxt = closest;
            auto node = 0;

            while(true) {
                if(closest < mint) {
                    break;
                }

                const auto& current = nodes[static_cast<size_t>(node)];
                if(current.leaf) {
                    auto limit = std::min(maxt, closest);
                    for(int i = 0; i < current.count; ++i) {
                        test(ordered[static_cast<size_t>(current.first + i)], mint, limit);
                        limit = std::min(maxt, closest);
                    }
                    if(stack.empty()) {
                        break;
                    }
                    const auto next = stack.back();
                    stack.pop_back();
                    node = next.node;
                    mint = next.mint;
                    maxt = next.maxt;
                    continue;
                }

                const auto axis = current.axis;
                const auto split = current.split;
                const auto t_plane = (split - r.origin[axis]) * (1.0f / r.direction[axis]);
                const bool left_of_split = r.origin[axis] < split
                    || (r.origin[axis] == split && r.direction[axis] <= 0.0f);
                const auto first = left_of_split ? current.left : current.right;
                const auto second = left_of_split ? current.right : current.left;

                if(t_plane > maxt || t_plane <= 0.0f) {
                    node = first;
                } else if(t_plane < mint) {
                    node = second;
                } else {
                    stack.push_back(item{second, t_plane, maxt});
                    node = first;
                    maxt = t_plane;
                }
            }
        }

    private:
        struct node {
            bool leaf = false;
            int axis = 0;
            float split = 0.0f;
            int left = -1;
            int right = -1;
            int first = 0;
            int count = 0;
        };

        struct edge_event {
            float position = 0.0f;
            std::uint32_t index = 0;
            std::uint8_t type = 0;
            std::uint8_t axis = 0;
        };

        static constexpr float traversal_cost = 15.0f;
        static constexpr float query_cost = 20.0f;
        static constexpr float empty_space_bonus = 0.9f;
        static constexpr int stop_prims = 6;
        static constexpr int max_bad_refines = 3;

        std::pair<global_precision, global_precision> split_probability(const aabb& box, int axis,
                                                                        global_precision left_width,
                                                                        global_precision right_width) const {
            const auto edge = box.maximum - box.minimum;
            const auto denom = edge.x * edge.y + edge.y * edge.z + edge.x * edge.z;
            const auto temp = 1.0f / denom;
            global_precision face = 0.0f;
            global_precision perimeter = 0.0f;
            if(axis == 0) {
                face = edge.y * edge.z;
                perimeter = edge.y + edge.z;
            } else if(axis == 1) {
                face = edge.x * edge.z;
                perimeter = edge.x + edge.z;
            } else {
                face = edge.x * edge.y;
                perimeter = edge.x + edge.y;
            }
            return {(face + perimeter * left_width) * temp, (face + perimeter * right_width) * temp};
        }

        void make_leaf(int index, const std::vector<std::uint32_t>& prims) {
            nodes[static_cast<size_t>(index)].leaf = true;
            nodes[static_cast<size_t>(index)].first = static_cast<int>(ordered.size());
            nodes[static_cast<size_t>(index)].count = static_cast<int>(prims.size());
            ordered.insert(ordered.end(), prims.begin(), prims.end());
        }

        float build_node(const std::vector<std::uint32_t>& prims, const aabb& voxel, int depth, int bad_refines) {
            const auto index = static_cast<int>(nodes.size());
            nodes.emplace_back();
            const auto prim_count = static_cast<int>(prims.size());
            const auto leaf_cost = static_cast<float>(prim_count) * query_cost;
            if(prim_count <= stop_prims || depth >= max_depth || voxel.surface_area() <= 0.0f) {
                make_leaf(index, prims);
                return leaf_cost;
            }

            std::vector<edge_event> events;
            events.reserve(static_cast<size_t>(prim_count) * 6);
            for(const auto prim : prims) {
                const auto& box = bounds[prim];
                for(std::uint8_t axis = 0; axis < 3; ++axis) {
                    const auto min_plane = box.minimum[axis];
                    const auto max_plane = box.maximum[axis];
                    if(min_plane == max_plane) {
                        events.push_back(edge_event{min_plane, prim, 1, axis});
                    } else {
                        events.push_back(edge_event{min_plane, prim, 2, axis});
                        events.push_back(edge_event{max_plane, prim, 0, axis});
                    }
                }
            }
            std::sort(events.begin(), events.end(), [](const edge_event& a, const edge_event& b) {
                if(a.axis != b.axis) {
                    return a.axis < b.axis;
                }
                if(a.position != b.position) {
                    return a.position < b.position;
                }
                if(a.type != b.type) {
                    return a.type < b.type;
                }
                return a.index < b.index;
            });

            int num_left[3] = {};
            int num_right[3] = {prim_count, prim_count, prim_count};
            auto best_cost = std::numeric_limits<float>::infinity();
            auto best_axis = 0;
            auto best_position = 0.0f;
            auto best_left = 0;
            auto best_right = 0;
            auto planar_left = false;

            for(size_t cursor = 0; cursor < events.size();) {
                const int axis = events[cursor].axis;
                const auto position = events[cursor].position;
                auto num_start = 0;
                auto num_end = 0;
                auto num_planar = 0;
                while(cursor < events.size() && events[cursor].axis == axis && events[cursor].position == position && events[cursor].type == 0) {
                    ++num_end;
                    ++cursor;
                }
                while(cursor < events.size() && events[cursor].axis == axis && events[cursor].position == position && events[cursor].type == 1) {
                    ++num_planar;
                    ++cursor;
                }
                while(cursor < events.size() && events[cursor].axis == axis && events[cursor].position == position && events[cursor].type == 2) {
                    ++num_start;
                    ++cursor;
                }

                num_right[axis] -= num_planar + num_end;
                if(position > voxel.minimum[axis] && position < voxel.maximum[axis]) {
                    const auto probability = split_probability(voxel, axis, position - voxel.minimum[axis], voxel.maximum[axis] - position);
                    const auto n_l = num_left[axis];
                    const auto n_r = num_right[axis];
                    if(num_planar == 0) {
                        auto cost = traversal_cost + query_cost * (probability.first * static_cast<float>(n_l) + probability.second * static_cast<float>(n_r));
                        if(n_l == 0 || n_r == 0) {
                            cost *= empty_space_bonus;
                        }
                        if(cost < best_cost) {
                            best_cost = cost;
                            best_axis = axis;
                            best_position = position;
                            best_left = n_l;
                            best_right = n_r;
                        }
                    } else {
                        auto cost_left = traversal_cost + query_cost * (probability.first * static_cast<float>(n_l + num_planar) + probability.second * static_cast<float>(n_r));
                        auto cost_right = traversal_cost + query_cost * (probability.first * static_cast<float>(n_l) + probability.second * static_cast<float>(n_r + num_planar));
                        if(n_l + num_planar == 0 || n_r == 0) {
                            cost_left *= empty_space_bonus;
                        }
                        if(n_l == 0 || n_r + num_planar == 0) {
                            cost_right *= empty_space_bonus;
                        }
                        if(cost_left < best_cost || cost_right < best_cost) {
                            planar_left = cost_left < cost_right;
                            best_cost = planar_left ? cost_left : cost_right;
                            best_axis = axis;
                            best_position = position;
                            best_left = planar_left ? n_l + num_planar : n_l;
                            best_right = planar_left ? n_r : n_r + num_planar;
                        }
                    }
                }
                num_left[axis] += num_start + num_planar;
            }

            if(best_cost >= leaf_cost) {
                if((best_cost > 4.0f * leaf_cost && prim_count < 16)
                    || bad_refines >= max_bad_refines
                    || best_cost == std::numeric_limits<float>::infinity()) {
                    make_leaf(index, prims);
                    return leaf_cost;
                }
                ++bad_refines;
            }

            std::vector<std::uint8_t> side(bounds.size(), 0);
            for(const auto prim : prims) {
                side[prim] = 0;
            }
            for(const auto& event : events) {
                if(event.axis != best_axis) {
                    continue;
                }
                if(event.type == 0 && event.position <= best_position) {
                    side[event.index] = 1;
                } else if(event.type == 2 && event.position >= best_position) {
                    side[event.index] = 2;
                } else if(event.type == 1) {
                    if(event.position < best_position || (event.position == best_position && planar_left)) {
                        side[event.index] = 1;
                    } else {
                        side[event.index] = 2;
                    }
                }
            }

            std::vector<std::uint32_t> left_prims;
            std::vector<std::uint32_t> right_prims;
            for(const auto prim : prims) {
                if(side[prim] == 1) {
                    left_prims.push_back(prim);
                } else if(side[prim] == 2) {
                    right_prims.push_back(prim);
                } else {
                    left_prims.push_back(prim);
                    right_prims.push_back(prim);
                }
            }

            aabb left_voxel = voxel;
            aabb right_voxel = voxel;
            left_voxel.maximum[best_axis] = best_position;
            right_voxel.minimum[best_axis] = best_position;

            const auto left_index = static_cast<int>(nodes.size());
            const auto left_cost = build_node(left_prims, left_voxel, depth + 1, bad_refines);
            const auto right_index = static_cast<int>(nodes.size());
            const auto right_cost = build_node(right_prims, right_voxel, depth + 1, bad_refines);
            const auto probability = split_probability(voxel, best_axis,
                best_position - voxel.minimum[best_axis], voxel.maximum[best_axis] - best_position);
            const auto final_cost = traversal_cost + probability.first * left_cost + probability.second * right_cost;
            if(final_cost < leaf_cost) {
                nodes[static_cast<size_t>(index)].leaf = false;
                nodes[static_cast<size_t>(index)].axis = best_axis;
                nodes[static_cast<size_t>(index)].split = best_position;
                nodes[static_cast<size_t>(index)].left = left_index;
                nodes[static_cast<size_t>(index)].right = right_index;
                return final_cost;
            }

            make_leaf(index, prims);
            return leaf_cost;
        }

        std::vector<aabb> bounds;
        std::vector<node> nodes;
        std::vector<std::uint32_t> ordered;
        int max_depth = 1;
    };
}
