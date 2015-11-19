#include "MVG.h"

//#include <ceres/ceres.h>

using namespace std;
using namespace Eigen;

Vector3d LinearLSTriangulation(const vector<Vector2d> &xs, const vector<pair<Matrix3d, Vector3d>> &states) {
    MatrixX3d A;
    VectorXd b;
    A.resize(xs.size() * 2, 3);
    b.resize(xs.size() * 2);
    size_t sstart = states.size() - xs.size();
    for (size_t i = 0; i < xs.size(); ++i) {
        const Vector2d & x = xs[i];
        const Matrix3d & R = states[sstart+i].first;
        const Vector3d & T = states[sstart+i].second;
        A.row(i * 2) = R.row(0) - x(0)*R.row(2);
        A.row(i * 2 + 1) = R.row(1) - x(1)*R.row(2);
        b(i * 2) = x(0)*T(2) - T(0);
        b(i * 2 + 1) = x(1)*T(2) - T(1);
    }
    
    return A.colPivHouseholderQr().solve(b);
}

Eigen::Vector3d LinearLSTriangulation(const std::vector<std::pair<Eigen::Vector2d, size_t>> &xs, const std::vector<std::pair<Eigen::Matrix3d, Eigen::Vector3d>> &states) {
    MatrixX3d A;
    VectorXd b;
    A.resize(xs.size() * 2, 3);
    b.resize(xs.size() * 2);
    for (size_t i = 0; i < xs.size(); ++i) {
        const Vector2d & x = xs[i].first;
        const Matrix3d & R = states[xs[i].second].first;
        const Vector3d & T = states[xs[i].second].second;
        A.row(i * 2) = R.row(0) - x(0)*R.row(2);
        A.row(i * 2 + 1) = R.row(1) - x(1)*R.row(2);
        b(i * 2) = x(0)*T(2) - T(0);
        b(i * 2 + 1) = x(1)*T(2) - T(1);
    }

    return A.colPivHouseholderQr().solve(b);
}

//Vector3d LinearLSTriangulation(const vector<Vector2d> &xs, const vector<Matrix3d> &Rs, const vector<Vector3d> &Ts) {
//    MatrixXd A;
//    VectorXd b;
//    A.resize(xs.size() * 2, 3);
//    b.resize(xs.size() * 2);
//    size_t u = xs.size();
//    size_t v = Rs.size();
//    size_t w = Ts.size();
//    size_t extent = min(u, min(v, w));
//    for (size_t i = 0; i < extent; ++i) {
//        u--; v--; w--;
//        auto & x = xs[u];
//        auto & R = Rs[v];
//        auto & T = Ts[w];
//        A.row(i * 2) = R.row(0) - x(0)*R.row(2);
//        A.row(i * 2 + 1) = R.row(1) - x(1)*R.row(2);
//        b(i * 2) = x(0)*T(2) - T(0);
//        b(i * 2 + 1) = x(1)*T(2) - T(1);
//    }
//
//    return A.colPivHouseholderQr().solve(b);
//}


//Vector3d LinearLSTriangulation(const vector<Vector2d> &xs, const vector<Matrix3d> &Rs, const vector<Vector3d> &Ts) {
//    MatrixXd A;
//    VectorXd b;
//    A.resize(xs.size() * 2, 3);
//    b.resize(xs.size() * 2);
//    size_t u = xs.size();
//    size_t v = Rs.size();
//    size_t w = Ts.size();
//    size_t extent = min(u, min(v, w));
//    for (size_t i = 0; i < extent; ++i) {
//        u--; v--; w--;
//        auto & x = xs[u];
//        auto & R = Rs[v];
//        auto & T = Ts[w];
//        A.row(i * 2) = R.row(0) - x(0)*R.row(2);
//        A.row(i * 2 + 1) = R.row(1) - x(1)*R.row(2);
//        b(i * 2) = x(0)*T(2) - T(0);
//        b(i * 2 + 1) = x(1)*T(2) - T(1);
//    }
//
//    return A.colPivHouseholderQr().solve(b);
//}
//
//struct TriangulateCostFunctor {
//    TriangulateCostFunctor(const Vector2d &x, const Matrix3d &r, const Vector3d &t) : x(x), r(r), t(t) {}
//
//    template <typename T>
//    bool operator() (const T* const p, T* residual) const {
//        T q[3];
//        for (int j = 0; j < 3; ++j) {
//            q[j] = T(t(j));
//            for (int i = 0; i < 3; ++i) {
//                q[j] += T(r(j, i))*p[i];
//            }
//        }
//        residual[0] = q[0] / q[2] - T(x[0]);
//        residual[1] = q[1] / q[2] - T(x[1]);
//        return true;
//    }
//
//    const Vector2d &x;
//    const Matrix3d &r;
//    const Vector3d &t;
//};
//
//Vector3d RefineTriangulation(const Vector3d &p0, const vector<Vector2d> &xs, const vector<Matrix3d> &Rs, const vector<Vector3d> &Ts) {
//    Vector3d p = p0;
//    ceres::Problem problem;
//    size_t u = xs.size();
//    size_t v = Rs.size();
//    size_t w = Ts.size();
//    size_t extent = min(u, min(v, w));
//    for (size_t i = 0; i < extent; ++i) {
//        u--; v--; w--;
//        ceres::CostFunction* cf = new ceres::AutoDiffCostFunction<TriangulateCostFunctor, 2, 3>(new TriangulateCostFunctor(xs[u], Rs[v], Ts[w]));
//        problem.AddResidualBlock(cf, nullptr, p.data());
//    }
//    ceres::Solver::Options options;
//    options.linear_solver_type = ceres::DENSE_NORMAL_CHOLESKY;
//    options.minimizer_progress_to_stdout = false;
//    ceres::Solver::Summary summary;
//    ceres::Solve(options, &problem, &summary);
//    return p;
//}
//
//Vector3d Triangulation(const vector<Vector2d> &xs, const vector<Matrix3d> &Rs, const vector<Vector3d> &Ts) {
//    Vector3d p0 = LinearLSTriangulation(xs, Rs, Ts);
//    return RefineTriangulation(p0, xs, Rs, Ts);
//}
