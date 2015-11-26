#include "MVG.h"

using namespace std;
using namespace Eigen;

Vector3d LinearLSTriangulation(const vector<Vector2d> &xs, const vector<CameraState> &states) {
    MatrixX3d A;
    VectorXd b;
    A.resize(xs.size() * 2, 3);
    b.resize(xs.size() * 2);
    size_t sstart = states.size() - xs.size();
    for (size_t i = 0; i < xs.size(); ++i) {
        const Vector2d & x = xs[i];
        const Matrix3d & R = states[sstart+i].R;
        const Vector3d & T = states[sstart+i].T;
        A.row(i * 2) = R.row(0) - x(0)*R.row(2);
        A.row(i * 2 + 1) = R.row(1) - x(1)*R.row(2);
        b(i * 2) = x(0)*T(2) - T(0);
        b(i * 2 + 1) = x(1)*T(2) - T(1);
    }
    
    return A.colPivHouseholderQr().solve(b);
}

Eigen::Vector3d LinearLSTriangulation(const vector<pair<Vector2d, size_t>> &xs, const vector<CameraState> &states) {
    MatrixX3d A;
    VectorXd b;
    A.resize(xs.size() * 2, 3);
    b.resize(xs.size() * 2);
    for (size_t i = 0; i < xs.size(); ++i) {
        const Vector2d & x = xs[i].first;
        const Matrix3d & R = states[xs[i].second].R;
        const Vector3d & T = states[xs[i].second].T;
        A.row(i * 2) = R.row(0) - x(0)*R.row(2);
        A.row(i * 2 + 1) = R.row(1) - x(1)*R.row(2);
        b(i * 2) = x(0)*T(2) - T(0);
        b(i * 2 + 1) = x(1)*T(2) - T(1);
    }

    return A.colPivHouseholderQr().solve(b);
}

Vector3d RefineTriangulation(const Vector3d &p0, const vector<pair<Vector2d, size_t>> &xs, const vector<CameraState> &states) {
    Vector3d p = p0;
    for (int iter = 0; iter < 10; ++iter) {
        Vector3d Jf = Vector3d::Zero();
        Matrix3d Hr = Matrix3d::Zero();
        for (size_t i = 0; i < xs.size(); ++i) {
            const Vector2d &x = xs[i].first;
            size_t ii = xs[i].second;
            const Matrix3d &R = states[ii].R;
            const Vector3d &T = states[ii].T;
            Vector3d pc = R*p + T;
            Vector2d xc(pc.x() / pc.z(), pc.y() / pc.z());
            MatrixXd Jpi(2, 3);
            Jpi.setIdentity();
            Jpi.col(2) = -xc;
            Jpi /= pc.z();
            Vector2d rc = xc - x;
            MatrixXd Jr = Jpi*R;
            Jf += rc.transpose()*Jr;
            Hr += Jr.transpose()*Jr;
        }
        Vector3d dp = Hr.colPivHouseholderQr().solve(Jf);
        p -= dp;
        if (dp.norm() < 1e-10) break;
    }

    return p;
}
