/*
    MVG.h - Multi-View Geometry
    多视图几何相关的函数
*/

#pragma once
#include <vector>
#include <Eigen/Eigen>
#include "JPL.h"


// 使用 Linear LS 方法进行三角化，参考[5]
// 该函数求解 p 使得下面的关系对于全部的 i 成立
//     Proj(Rs[i]*p+Ts[i]) = xs[i]
// xs 和 states 从后向前对应，xs 的长度不超过 states 的长度
Eigen::Vector3d LinearLSTriangulation(const std::vector<Eigen::Vector2d> &xs, const std::vector<std::pair<Eigen::Matrix3d, Eigen::Vector3d>> &states);

// 使用 Linear LS 方法进行三角化，参考[5]
// 该函数求解 p 使得下面的关系对于全部的 i 成立
//     Proj(Rs[i]*p+Ts[i]) = xs[i]
// xs.second 代表 states 中的对应项
Eigen::Vector3d LinearLSTriangulation(const std::vector<std::pair<Eigen::Vector2d, size_t>> &xs, const std::vector<std::pair<Eigen::Matrix3d, Eigen::Vector3d>> &states);

// 根据初始值 p0 对三角化进行最小二乘优化
Eigen::Vector3d RefineTriangulation(const Eigen::Vector3d &p0, const std::vector<std::pair<Eigen::Vector2d, size_t>> &xs, const std::vector<std::pair<Eigen::Matrix3d, Eigen::Vector3d>> &states);
