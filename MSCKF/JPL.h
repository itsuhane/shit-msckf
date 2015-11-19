/*
    JPL.h - JPL's Quaternion
    各种适配 JPL 四元数运算的函数，参考 [3]
*/

#pragma once
#include <Eigen/Eigen>

// JPL 四元数实际上使用一个四维向量表示
typedef Eigen::Vector4d JPL_Quaternion;

// 标准化一个 JPL 四元数
// 为了表示旋转，四元数需要为单位长
// 同时为了更好的数值稳定性，保证 w 大于零
inline JPL_Quaternion JPL_Normalize(const JPL_Quaternion &q) {
    JPL_Quaternion r = q.normalized();
    if (r.w() < 0) r = -r;
    return r;
}

// 将一个 JPL 四元数转换成 Hamilton 四元数
inline Eigen::Quaterniond JPL_toHamilton(const JPL_Quaternion &q) {
    return Eigen::Quaterniond(q.w(), -q.x(), -q.y(), -q.z());
}

// 将一个 Hamilton 四元数转换为 JPL 四元数
inline JPL_Quaternion HamiltonToJPL(const Eigen::Quaterniond& q) {
    return JPL_Quaternion(-q.x(), -q.y(), -q.z(), q.w());
}

// 获得一个 JPL 四元数对应的旋转
inline Eigen::Matrix3d JPL_C(const JPL_Quaternion &q) {
    return JPL_toHamilton(q).toRotationMatrix();
}

// 获得一个 JPL 四元数对应的旋转的转置
inline Eigen::Matrix3d JPL_CT(const JPL_Quaternion &q) {
    // 先计算共轭，再求旋转矩阵
    return Eigen::Quaterniond(q.w(), q.x(), q.y(), q.z()).toRotationMatrix();
}

// 计算两个 JPL 四元数的乘积
inline JPL_Quaternion JPL_Multiply(const JPL_Quaternion &q, const JPL_Quaternion &p) {
    return JPL_Quaternion(
        q(3)*p(0) + q(2)*p(1) - q(1)*p(2) + q(0)*p(3),
        -q(2)*p(0) + q(3)*p(1) + q(0)*p(2) + q(1)*p(3),
        q(1)*p(0) - q(0)*p(1) + q(3)*p(2) + q(2)*p(3),
        -q(0)*p(0) - q(1)*p(1) - q(2)*p(2) + q(3)*p(3)
        );
}

// [3] 中定义的 [ * x] 算子
inline Eigen::Matrix3d JPL_Cross(const Eigen::Vector3d &w) {
    return (Eigen::Matrix3d() << 0.0, -w.z(), w.y(),
        w.z(), 0.0, -w.x(),
        -w.y(), w.x(), 0.0).finished();
}

// [3] 中定义的 Omega 算子
inline Eigen::Matrix4d JPL_Omega(const Eigen::Vector3d& w) {
    return (Eigen::Matrix4d() << 0.0, w.z(), -w.y(), w.x(),
        -w.z(), 0.0, w.x(), w.y(),
        w.y(), -w.x(), 0.0, w.z(),
        -w.x(), -w.y(), -w.z(), 0.0).finished();
}