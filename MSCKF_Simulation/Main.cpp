#include <iostream>
#include <functional>
#include "MSCKF.h"
#include <qb/meshio.hpp>

using namespace std;
using namespace Eigen;

void Simulation_Propagation_Static(); // 静止不动
void Simulation_Propagation_FixedSpeed(); // 匀速直线运动
void Simulation_Propagation_FixedAcceleration(); // 匀加速直线
void Simulation_Propagation_SinusoidalShifting(); // 2D正弦曲线振荡
void Simulation_Propagation_CircularShifting(); // 圆形平移
void Simulation_Propagation_CircularDriving(); // 圆形前进

int main(int argc, char* argv[]) {
    Simulation_Propagation_Static();
    Simulation_Propagation_FixedSpeed();
    Simulation_Propagation_FixedAcceleration();
    Simulation_Propagation_SinusoidalShifting();
    Simulation_Propagation_CircularShifting();
    Simulation_Propagation_CircularDriving();
    return 0;
}

void addAxis(qb::mesh &m, const Quaterniond &orientation, const Vector3d &position, bool highlight = true) {
    Matrix3d O = orientation.toRotationMatrix();
    Vector3d X = O.col(0);
    Vector3d Y = O.col(1);
    Vector3d Z = O.col(2);
    unsigned char v = highlight ? 255 : 192;
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
