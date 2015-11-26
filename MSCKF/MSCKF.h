/*
    MSCKF.h - Multi-State Constrained Kalman Filter
    核心的滤波器模块
*/

#pragma once

#include <vector>
#include <unordered_map>
#include <functional>
#include <Eigen/Eigen>
#include "JPL.h"
#include "MVG.h"

namespace Eigen {
    typedef Matrix<double, 15, 15> Matrix15d; // 15x15 矩阵
}

class MSCKF {
public:
    MSCKF();

    // 设置噪音的协方差大小
    void setNoiseCov(const Eigen::Matrix3d &cov_ng, const Eigen::Matrix3d &cov_nwg, const Eigen::Matrix3d &cov_na, const Eigen::Matrix3d &cov_nwa, double sigma_im_squared);

    // 初始化滤波器
    // 需要提供初始的旋转、角速度偏差、速度、加速度偏差、位置、重力加速度
    void initialize(const JPL_Quaternion &q, const Eigen::Vector3d &bg, const Eigen::Vector3d &v, const Eigen::Vector3d &ba, const Eigen::Vector3d &p, double g = 9.81);

    // 用 t 时刻的角速度和加速度来进行积分
    void propagate(double t, const Eigen::Vector3d &w, const Eigen::Vector3d &a);

    // 使用 GPS 信号进行更新
    void track(double t, const Eigen::Vector3d &position, double sigma_p);

    // 进行特征点跟踪，并进行增广或矫正
    // features[跟踪到的新帧中的特征id] = <被跟踪到的老帧中的特征id，特征在新帧中的投影位置（已经除掉内参）>
    void track(double t, const std::unordered_map<size_t, std::pair<size_t, Eigen::Vector2d>> &matches);

    // 获得当前传感器在世界坐标系中的朝向
    Eigen::Quaterniond orientation() const { return JPL_toHamilton(m_q).conjugate(); }

    // 获得当前传感器在世界坐标系中的位置
    Eigen::Vector3d position() const { return m_p; }

    Eigen::Matrix3d positionCovariance() const { return m_PII.block<3, 3>(12, 12); }

    // 获得当前相机在世界坐标系中的朝向
    Eigen::Quaterniond cameraOrientation() const { return JPL_toHamilton(JPL_Multiply(m_q_imu_to_cam, m_q)).conjugate(); }

    // 获得当前相机在世界坐标系中的位置
    Eigen::Vector3d cameraPosition() const { return m_p + JPL_CT(m_q)*m_p_cam_in_imu; }

private:
    // MotionSystem 需要在积分时获得当前的各种参数
    friend class MotionSystem;

    size_t m_state_limit;            // 相机状态的个数上限
    JPL_Quaternion m_q_imu_to_cam;   // 从IMU坐标系到相机坐标系的旋转
    Eigen::Vector3d m_p_cam_in_imu;  // 相机中心在IMU坐标系中的坐标
    double m_sigma_im_squared;       // 投影坐标的均方差

    Eigen::Matrix3d m_cov_ng;        // 角速度传感器协方差
    Eigen::Matrix3d m_cov_nwg;       // 角速度偏移噪音协方差
    Eigen::Matrix3d m_cov_na;        // 加速度传感器协方差
    Eigen::Matrix3d m_cov_nwa;       // 加速度偏移噪音协方差

    JPL_Quaternion m_q;              // 世界坐标系到设备坐标系的旋转
    Eigen::Vector3d m_bg;            // 设备坐标系中角速度偏差
    Eigen::Vector3d m_v;             // 世界坐标系中的速度
    Eigen::Vector3d m_ba;            // 设备坐标系中加速度偏差
    Eigen::Vector3d m_p;             // 世界坐标系中的位置
    double m_g;                      // 重力加速度大小

    Eigen::Matrix15d m_PII;          // IMU状态的协方差
    std::vector<Eigen::MatrixXd> m_PIC; // 分块表示的 PIC (IMU和相机状态之间的协方差)
    std::vector<std::vector<Eigen::MatrixXd>> m_PCC; // 分块表示的 PCC 的上三角 (相机状态的协方差)

    bool m_has_old = false;
    double m_t_old;
    Eigen::Vector3d m_w_old;
    Eigen::Vector3d m_a_old;

    std::vector<CameraState> m_states; // 相机 state，我们使用 R 和 T 代替论文中的 q 和 p
    std::unordered_map<size_t, std::vector<Eigen::Vector2d>> m_tracks; // 特征 tracks
};
