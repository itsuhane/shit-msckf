#pragma once
#include <vector>
#include <Eigen/Eigen>

class GyroCalibrator {
public:
    void sensor(double t, const Eigen::Vector3d &w, const Eigen::Vector3d &a);
    bool vibrated() const;
    bool calibrated() const { return m_calibrated; }
    Eigen::Vector3d bw() const { return m_bw; }
    void reset();

private:
    std::vector<Eigen::Vector3d> m_ws;
    std::vector<Eigen::Vector3d> m_as;

    bool m_calibrated = false;
    Eigen::Vector3d m_bw = Eigen::Vector3d::Zero();
};
