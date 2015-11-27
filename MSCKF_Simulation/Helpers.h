#pragma once

#include <random>
#include <Eigen/Eigen>
#include <qb/meshio.hpp>

struct Uniform {
    std::mt19937 gen;
    std::uniform_real_distribution<> d;

    Uniform(double minimum = 0.0, double maximum = 1.0, int seed = 0) : gen(seed), d(minimum, maximum) {}

    double value() {
        return d(gen);
    }

    Eigen::Vector2d vector2() {
        return Eigen::Vector2d(value(), value());
    }

    Eigen::Vector3d vector3() {
        return Eigen::Vector3d(value(), value(), value());
    }
};

struct Gaussian {
    std::mt19937 gen;
    std::normal_distribution<> d;

    Gaussian(double sigma = 1.0, int seed = 0) : gen(seed), d(0, sigma) {}

    double value() {
        return d(gen);
    }

    Eigen::Vector2d vector2() {
        return Eigen::Vector2d(value(), value());
    }

    Eigen::Vector3d vector3() {
        return Eigen::Vector3d(value(), value(), value());
    }
};

void addLine(qb::mesh &m, const Eigen::Vector3d &p, const Eigen::Vector3d &d, const std::array<unsigned char, 3> &c) {
    for (int i = 0; i < 100; ++i) {
        m.vertices.push_back(p + d*i / 100.0);
        m.colors.push_back(c);
    }
}

void addAxis(qb::mesh &m, const Eigen::Quaterniond &orientation, const Eigen::Vector3d &position, unsigned char brightness = 255) {
    Eigen::Matrix3d O = orientation.toRotationMatrix();
    Eigen::Vector3d X = O.col(0);
    Eigen::Vector3d Y = O.col(1);
    Eigen::Vector3d Z = O.col(2);
    addLine(m, position, 0.1*X, { brightness, 0, 0 });
    addLine(m, position, 0.1*Y, { 0, brightness, 0 });
    addLine(m, position, 0.1*Z, { 0, 0, brightness });
}
