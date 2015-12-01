#pragma once

#include <vector>
#include <unordered_map>
#include <memory>
#include <Eigen/Eigen>
#include <SiftGPU/SiftGPU.h>

class SiftTracker {
public:
    SiftTracker();

    std::unordered_map<size_t, std::pair<size_t, Eigen::Vector2d>> track(std::vector<unsigned char> &img);

private:
    int m_img_width;
    int m_img_height;

    double m_camera_fx;
    double m_camera_fy;
    double m_camera_cx;
    double m_camera_cy;

    std::unique_ptr<SiftGPU> m_sift;
    std::unique_ptr<SiftMatchGPU> m_matcher;
};
