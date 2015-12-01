#include "AccelerometerCalibrator.h"
#include <iostream>

using namespace std;
using namespace Eigen;

AccelerometerCalibrator::AccelerometerCalibrator() {
    m_references.emplace_back(0, 0, 1);
    m_references.emplace_back(0, 0, -1);
    m_references.emplace_back(-1, 0, 0);
    m_references.emplace_back(0, -1, 0);
    m_references.emplace_back(1, 0, 0);
    m_references.emplace_back(0, 1, 0);
}

void AccelerometerCalibrator::sensor(double t, const Vector3d &w, const Vector3d &a) {
    if (m_calibrated) {
        return;
    }

    m_ws.push_back(w);
    m_as.push_back(a);

    const Vector3d &current_reference = m_references[m_rps.size()];
    if (vibrated() || m_as.back().dot(current_reference)<0.9*m_as.back().norm()) {
        resetrp();
    }

    cout << m_as.size() << endl;

    if (m_as.size() == 100) {
        Vector3d rp = Vector3d::Zero();
        for (auto & a : m_as) {
            rp += a;
        }
        rp /= (double)m_as.size();
        m_rps.push_back(rp);
        cout << "Anchor" << endl;
    }

    if (m_rps.size() == m_references.size()) {

        Vector3d ba = Vector3d::Zero();
        double g = 9.81;

        for (int iter = 0; iter < 10; ++iter) {
            Matrix4d A = Matrix4d::Zero();
            Vector4d b = Vector4d::Zero();
            for (size_t i = 0; i < m_rps.size(); ++i) {
                double rg = (m_rps[i] - ba).norm();
                double ri = rg - g;
                Vector3d dridba = Vector3d::Zero();
                if (rg >= 1e-7) {
                    dridba = (ba - m_rps[i]) / rg;
                }
                Vector4d dri(dridba.x(), dridba.y(), dridba.z(), -1);
                Matrix4d Ari = dri*dri.transpose();
                Vector4d bri = dri*ri;
                A += Ari;
                b += bri;
            }
            Vector4d dx = A.colPivHouseholderQr().solve(b);
            ba.x() -= dx.x();
            ba.y() -= dx.y();
            ba.z() -= dx.z();
            g -= dx.w();
            if (dx.norm() < 1e-7) break;
        }
        cout << "ba = " << ba.transpose() << endl;
        cout << "g = " << g << endl;
        m_calibrated = true;
    }
}

bool AccelerometerCalibrator::vibrated() const {
    if (m_as.size() > 1) {
        return (m_as[m_as.size() - 1] - m_as[m_as.size() - 2]).norm() > 0.7;
    }
    else {
        return false;
    }
}

//void AccelerometerCalibrator::reset();

void AccelerometerCalibrator::resetrp() {
    m_ws.clear();
    m_as.clear();
}