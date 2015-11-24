#include <iostream>
#include <functional>
#include <random>
#include "MSCKF.h"
#include <qb/meshio.hpp>

using namespace std;
using namespace Eigen;

struct Gaussian {
    random_device rd;
    mt19937 gen;
    normal_distribution<> d;

    Gaussian(double sigma) : rd(), gen(rd()), d(0, sigma) {}

    double next() {
        return d(gen);
    }

    Vector3d vector() {
        return Vector3d(next(), next(), next());
    }
};

void Simulation_Propagation_Static(); // 静止不动
void Simulation_Propagation_FixedSpeed(); // 匀速直线运动
void Simulation_Propagation_FixedAcceleration(); // 匀加速直线
void Simulation_Propagation_SinusoidalShifting(); // 2D正弦曲线振荡
void Simulation_Propagation_CircularShifting(); // 圆形平移
void Simulation_Propagation_CircularDriving(); // 圆形前进

void Simulation_CircularDrivingTest();

int main(int argc, char* argv[]) {
    Simulation_CircularDrivingTest();
    return 0;
}

void addLine(qb::mesh &m, const Vector3d &p, const Vector3d &d) {
    for (int i = 0; i < 100; ++i) {
        m.vertices.push_back(p+d*i/99.0);
        m.colors.push_back({ 255, 255, 0 });
    }
}

void addAxis(qb::mesh &m, const Quaterniond &orientation, const Vector3d &position, bool highlight = true) {
    Matrix3d O = orientation.toRotationMatrix();
    Vector3d X = O.col(0);
    Vector3d Y = O.col(1);
    Vector3d Z = O.col(2);
    unsigned char v = highlight ? 255 : 96;
    for (int i = 0; i < 100; ++i) {
        m.vertices.push_back(position + 0.001*X*i);
        m.colors.push_back({ v, 0, 0 });
        m.vertices.push_back(position + 0.001*Y*i);
        m.colors.push_back({ 0, v, 0 });
        m.vertices.push_back(position + 0.001*Z*i);
        m.colors.push_back({ 0, 0, v });
    }
}

void Simulate(MSCKF &ekf, qb::mesh &output, function<Vector3d(double)> gyro, function<Vector3d(double)> accelerometer) {
    for (double t = 0.0; t <= 10.0; t += 0.01) {
        ekf.propagate(t, gyro(t), accelerometer(t));
        addAxis(output, ekf.orientation(), ekf.position());
    }
}

void Simulation_Propagation_Static() {
    qb::mesh output;
    MSCKF ekf;
    ekf.initialize(JPL_Quaternion::UnitW(), Vector3d::Zero(), Vector3d::Zero(), Vector3d::Zero(), Vector3d::Zero(), 0.0);
    Simulate(ekf, output,
        [](double t) -> Vector3d { return Vector3d::Zero(); },
        [](double t) -> Vector3d { return Vector3d::Zero(); }
        );
    output.save("E:\\SensorFusion\\Propagation-Static-IMU.ply");
}

void Simulation_Propagation_FixedSpeed() {
    qb::mesh output;
    MSCKF ekf;
    ekf.initialize(JPL_Quaternion::UnitW(), Vector3d::Zero(), Vector3d(0.1, 0.0, 0.0), Vector3d::Zero(), Vector3d::Zero(), 0.0);
    Simulate(ekf, output,
        [](double t) -> Vector3d { return Vector3d::Zero(); },
        [](double t) -> Vector3d { return Vector3d::Zero(); }
    );
    output.save("E:\\SensorFusion\\Propagation-FixedSpeed-IMU.ply");
}

void Simulation_Propagation_FixedAcceleration() {
    qb::mesh output;
    MSCKF ekf;
    ekf.initialize(JPL_Quaternion::UnitW(), Vector3d::Zero(), Vector3d::Zero(), Vector3d::Zero(), Vector3d::Zero(), 0.0);
    Simulate(ekf, output,
        [](double t) -> Vector3d { return Vector3d::Zero(); },
        [](double t) -> Vector3d { return Vector3d(0.01, 0.0, 0.0); }
    );
    output.save("E:\\SensorFusion\\Propagation-FixedAcceleration-IMU.ply");
}

void Simulation_Propagation_SinusoidalShifting() {
    const double R = 1.0;
    const double T = 5.0;
    const double w = 2.0 * M_PI / T;
    qb::mesh output;
    MSCKF ekf;
    ekf.initialize(JPL_Quaternion::UnitW(), Vector3d::Zero(), Vector3d(R, R*w, 0.0), Vector3d::Zero(), Vector3d::Zero(), 0.0);
    Simulate(ekf, output,
        [](double t) -> Vector3d { return Vector3d::Zero(); },
        [=](double t) -> Vector3d { return Vector3d(0.0, -R*w*w*sin(w*t), 0.0); }
    );
    output.save("E:\\SensorFusion\\Propagation-SinusoidalShifting-IMU.ply");
}

void Simulation_Propagation_CircularShifting() {
    const double R = 1.0;
    const double T = 5.0;
    const double w = 2.0 * M_PI / T;
    qb::mesh output;
    MSCKF ekf;
    ekf.initialize(JPL_Quaternion::UnitW(), Vector3d::Zero(), Vector3d(R*w, 0.0, 0.0), Vector3d::Zero(), Vector3d::Zero(), 0.0);
    Simulate(ekf, output,
        [](double t) -> Vector3d { return Vector3d::Zero(); },
        [=](double t) -> Vector3d { return Vector3d(-R*w*w*sin(w*t), R*w*w*cos(w*t), 0.0); }
    );
    output.save("E:\\SensorFusion\\Propagation-CircularShifting-IMU.ply");
}

void Simulation_Propagation_CircularDriving() {
    const double R = 1.0;
    const double T = 5.0;
    const double w = 2.0 * M_PI / T;
    qb::mesh output;
    MSCKF ekf;
    ekf.initialize(JPL_Quaternion::UnitW(), Vector3d::Zero(), Vector3d(R*w, 0.0, 0.0), Vector3d::Zero(), Vector3d::Zero(), 0.0);
    Simulate(ekf, output,
        [=](double t) -> Vector3d { return Vector3d(0.0, 0.0, w); },
        [=](double t) -> Vector3d { return Vector3d(0.0, R*w*w, 0.0); }
    );
    output.save("E:\\SensorFusion\\Propagation-CircularDriving-IMU.ply");
}

void Simulation_CircularDrivingTest() {
    const double R = 5.0;
    const double T = 5.0;
    const double w = 2.0 * M_PI / T;
    
    Quaterniond InitOrientation((Matrix3d() << Vector3d::UnitX(), -Vector3d::UnitZ(), Vector3d::UnitY()).finished());

    qb::mesh features;
    for (int i = 0; i < 360; i += 5) {
        for (int j = -10; j < 10; j+=2) {
            double angle = i*M_PI / 180.0;
            double height = j*0.1;
            Vector3d point(4.0*cos(angle), 4.0*sin(angle), height);
            features.vertices.push_back(point);
        }
    }
    //features.save("E:\\SensorFusion\\Synthesis\\CircularDrivingTest-Features.ply");

    qb::mesh raw, error, truth;
    MSCKF rkf, ekf, tkf;
    rkf.initialize(HamiltonToJPL(InitOrientation), Vector3d::Zero(), Vector3d(R*w, 0.0, 0.0), Vector3d::Zero(), Vector3d(0, -R, 0), 0.0);
    ekf.initialize(HamiltonToJPL(InitOrientation), Vector3d::Zero(), Vector3d(R*w, 0.0, 0.0), Vector3d::Zero(), Vector3d(0, -R, 0), 0.0);
    tkf.initialize(HamiltonToJPL(InitOrientation), Vector3d::Zero(), Vector3d(R*w, 0.0, 0.0), Vector3d::Zero(), Vector3d(0, -R, 0), 0.0);

    int count = 0;
    for (double t = 0.0; t <= 10.0; t += 0.01) {
        Gaussian g(0.1);
        rkf.setNoiseCov(Matrix3d::Identity()*0.1, Matrix3d::Identity()*0.0, Matrix3d::Identity()*0.1, Matrix3d::Identity()*0.0, 0.0);
        ekf.setNoiseCov(Matrix3d::Identity()*0.1, Matrix3d::Identity()*0.0, Matrix3d::Identity()*0.1, Matrix3d::Identity()*0.0, 0.0);

        Vector3d gyro(0.0, w, 0.0);
        Vector3d acce(0.0, 0.0, -R*w*w);
        Vector3d gyro_wrong = gyro + g.vector();
        Vector3d acce_wrong = acce + g.vector();
        rkf.propagate(t, gyro_wrong, acce_wrong);
        ekf.propagate(t, gyro_wrong, acce_wrong);
        tkf.propagate(t, gyro, acce);
        addAxis(raw, rkf.orientation(), rkf.position(), false);
        addAxis(truth, tkf.orientation(), tkf.position(), false);
        if (count%4==0) {
#if 1
            qb::mesh frame;
            Matrix3d R = tkf.cameraOrientation().toRotationMatrix().transpose();
            Vector3d P = tkf.cameraPosition();
            unordered_map<size_t, pair<size_t, Vector2d>> matches;
            for (size_t i = 0; i < features.vertices.size(); ++i) {
                Vector3d p = features.vertices[i];
                Vector3d q = R*(p - P);
                if (q.z() < 2.0 && abs(q.x()) < 1.0 && abs(q.y()) < 1.0) {
                    frame.vertices.push_back(p);
                    matches[i] = make_pair(i, Vector2d(q.x()/q.z(), q.y()/q.z()));
                }
            }
            ekf.track(t, matches);
            if (ekf.position().norm() < 20) {
                addAxis(error, ekf.orientation(), ekf.position());
                if (ekf.dp.norm() < 5) {
                    addLine(error, ekf.position(), ekf.dp);
                }
            }
#endif
        }
        count++;
    }

    raw.save("E:\\SensorFusion\\Synthesis\\CircularDrivingTest-Raw.ply");
    truth.save("E:\\SensorFusion\\Synthesis\\CircularDrivingTest-Truth.ply");
    error.save("E:\\SensorFusion\\Synthesis\\CircularDrivingTest-Error.ply");
}
