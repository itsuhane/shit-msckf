#include "GyroCalibrator.h"

using namespace std;
using namespace Eigen;

void GyroCalibrator::sensor(double t, const Vector3d &w, const Vector3d &a) {
    if (m_calibrated) {
        return;
    }

    m_ws.push_back(w);
    m_as.push_back(a);
    if (vibrated()) {
        reset();
    }

    if (m_ws.size() == 1000) {
        m_bw = Vector3d::Zero();
        for (auto & w : m_ws) {
            m_bw += w;
        }
        m_bw /= (double)m_ws.size();
        m_calibrated = true;
    }
}

bool GyroCalibrator::vibrated() const {
    if (m_as.size() > 1) {
        return (m_as[m_as.size() - 1] - m_as[m_as.size() - 2]).norm() > 0.7;
    }
    else {
        return false;
    }
}

void GyroCalibrator::reset() {
    m_ws.clear();
    m_as.clear();
    m_calibrated = false;
    m_bw = Vector3d::Zero();
}
