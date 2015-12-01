#pragma once
#include <vector>
#include <Eigen/Eigen>

class AccelerometerCalibrator {
public:
    AccelerometerCalibrator();

    void sensor(double t, const Eigen::Vector3d &w, const Eigen::Vector3d &a);
    bool vibrated() const;
    bool calibrated() const { return m_calibrated; }
    Eigen::Vector3d ba() const { return m_ba; }
    void reset();

private:
    void resetrp();

    std::vector<Eigen::Vector3d> m_references;

    std::vector<Eigen::Vector3d> m_ws;
    std::vector<Eigen::Vector3d> m_as;
    std::vector<Eigen::Vector3d> m_rps;

    bool m_calibrated = false;
    Eigen::Vector3d m_ba = Eigen::Vector3d::Zero();
};