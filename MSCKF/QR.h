/*
    QR.h - QR Decomposition Related
    一些利用 QR 分解的工具，用来做不显式分解的子空间投影
*/

#pragma once
#include <Eigen/Eigen>

// 计算 Givens 旋转用到的系数
// 其中 a 为需要置零的元素， b 为用来置零的元素
// 返回 c 为对应的 cos 分量，s 为 sin 分量
inline void Givens(double a, double b, double &c, double &s) {
    if (abs(b)<1.0e-15) {
        c = copysign(1.0, a);
        s = 0;
    }
    else if (abs(a)<1.0e-15) {
        c = 0;
        s = -copysign(1.0, b);
    }
    else if (abs(b) > abs(a)) {
        double t = a / b;
        double u = copysign(sqrt(1 + t*t), b);
        s = -1.0 / u;
        c = -s*t;
    }
    else {
        double t = b / a;
        double u = copysign(sqrt(1 + t*t), a);
        c = 1.0 / u;
        s = -c*t;
    }
}

