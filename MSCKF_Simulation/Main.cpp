#include <fstream>
#include <iostream>
#include <functional>
#include <boost/lexical_cast.hpp>
#include "MSCKF.h"

#include "Helpers.h"

using namespace std;
using namespace Eigen;

const double R = 5.0;
Gaussian noise(1.0, 10);

vector<Vector3d> makeFeatures() {
    vector<Vector3d> features;
    Uniform rnd(-1.0, 1.0);
    for (int i = 0; i < 500; ++i) {
        Vector3d p = rnd.vector3()*R;
        if (abs(p.z()) < 3.0 && Vector2d(p.x(), p.y()).norm() < R) {
            features.push_back(p);
        }
    }
    return features;
}

void saveFeatures(const std::string &filepath, vector<Vector3d> features) {
    qb::mesh m;
    m.vertices.swap(features);
    m.save(filepath);
}

unordered_map<size_t, pair<size_t, Vector2d>> frameFeature(const vector<Vector3d> &features, const Quaterniond &orientation, const Vector3d &position, double nim) {
    unordered_map<size_t, pair<size_t, Vector2d>> frame;
    Matrix3d R = orientation.toRotationMatrix().transpose();
    for (size_t i = 0; i < features.size(); ++i) {
        const Vector3d &p = features[i];
        Vector3d q = R*(p - position);
        Vector2d z(q.x() / q.z(), q.y() / q.z());
        if (q.z()>0.5 && q.z() < 4.0 && abs(z.x()) < 1.0 && abs(z.y()) < 1.0) {
            frame[i] = make_pair(i, z+noise.vector2()*nim);
        }
    }
    return frame;
}

void saveFrame(const std::string &filepath, const vector<Vector3d> &features, const unordered_map<size_t, pair<size_t, Vector2d>> &frame) {
    qb::mesh m;
    for (auto &f : frame) {
        m.vertices.push_back(features[f.first]);
    }
    m.save(filepath);
}

int main(int argc, char* argv[]) {
    const double T = 25.0;
    const double w = 2.0 * M_PI / T;
    double ng = 0.002;
    double nwg = 0.0005;
    double na = 0.05;
    double nwa = 0.005;
    double nim = 0.005;

    //nim = boost::lexical_cast<double>(argv[1]);
    //cout << "nim = " << nim << endl;

    vector<Vector3d> features = makeFeatures();

    Quaterniond initOrientation((Matrix3d() << Vector3d::UnitX(), -Vector3d::UnitZ(), Vector3d::UnitY()).finished());
    Vector3d gyro(0.0, w, 0.0);
    Vector3d acce(0.0, 0.0, -R*w*w);

    MSCKF ekf, ekf_propagate_only, ekf_true;
    qb::mesh ekf_path, ekf_propagate_only_path, ekf_true_path;

    ekf.setNoiseCov(Matrix3d::Identity()*ng, Matrix3d::Identity()*nwg, Matrix3d::Identity()*na, Matrix3d::Identity()*nwa, nim);
    ekf_propagate_only.setNoiseCov(Matrix3d::Identity()*ng, Matrix3d::Identity()*nwg, Matrix3d::Identity()*na, Matrix3d::Identity()*nwa, nim);
    ekf_true.setNoiseCov(Matrix3d::Identity()*ng, Matrix3d::Identity()*nwg, Matrix3d::Identity()*na, Matrix3d::Identity()*nwa, nim);
    ekf.initialize(HamiltonToJPL(initOrientation), Vector3d::Zero(), Vector3d(R*w, 0.0, 0.0), Vector3d::Zero(), Vector3d(0, -R, 0), 0.0);
    ekf_propagate_only.initialize(HamiltonToJPL(initOrientation), Vector3d::Zero(), Vector3d(R*w, 0.0, 0.0), Vector3d::Zero(), Vector3d(0, -R, 0), 0.0);
    ekf_true.initialize(HamiltonToJPL(initOrientation), Vector3d::Zero(), Vector3d(R*w, 0.0, 0.0), Vector3d::Zero(), Vector3d(0, -R, 0), 0.0);

    ofstream log("E:/SensorFusion/Synthesis/log-" + boost::lexical_cast<string>(nim) + ".csv", ofstream::out | ofstream::trunc);

    for (double t = 0.0; t <= 50.0; t += 0.01) {
        Vector3d noisy_gyro = gyro + noise.vector3()*ng;
        Vector3d noisy_acce = acce + noise.vector3()*na + Vector3d(0.02, 0.2, 0.001);
        ekf.propagate(t, noisy_gyro, noisy_acce);
        ekf_propagate_only.propagate(t, noisy_gyro, noisy_acce);
        ekf_true.propagate(t, gyro, acce);
        size_t fsize = 0;
        auto frame = frameFeature(features, ekf_true.cameraOrientation(), ekf_true.cameraPosition(), nim);
        ekf.update(t, frame);
        fsize = frame.size();
        addAxis(ekf_path, ekf.orientation(), ekf.position());
        addAxis(ekf_propagate_only_path, ekf_propagate_only.orientation(), ekf_propagate_only.position());
        addAxis(ekf_true_path, ekf_true.orientation(), ekf_true.position());
        stringstream s;
        s << t << ", " << fsize << ", " << (ekf_propagate_only.position() - ekf_true.position()).norm() << ", " << (ekf.position() - ekf_true.position()).norm() << ", " << ekf.positionCovariance().diagonal().sum();
        cout << s.str() << endl;
        log << s.str() << endl;
    }

    ekf_path.save("E:/SensorFusion/Synthesis/path-" + boost::lexical_cast<string>(nim) + ".ply");
    ekf_propagate_only_path.save("E:/SensorFusion/Synthesis/path-propagate.ply");
    ekf_true_path.save("E:/SensorFusion/Synthesis/path-true.ply");

    return 0;
}
