#include "SiftTracker.h"

#include <array>
#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <Windows.h>
#include <gl/GL.h>

using namespace std;
using namespace Eigen;

SiftTracker::SiftTracker() {
    m_img_width = 640;
    m_img_height = 480;
    m_camera_fx = 613.969482;
    m_camera_fy = 613.969482;
    m_camera_cx = 320.0;
    m_camera_cy = 240.0;

    m_sift.reset(new SiftGPU);
    char * siftargv[] = { "-fo", "-1", "-v", "0", "-noprep" };
    m_sift->ParseParam(5, siftargv);
    if (m_sift->CreateContextGL() != SiftGPU::SIFTGPU_FULL_SUPPORTED) {
        printf("Cannot create SiftGPU.\n");
        throw "Cannot create SiftGPU.";
    }

    m_matcher.reset(new SiftMatchGPU(4096));
    if (m_matcher->VerifyContextGL() == 0) {
        printf("Cannot create SiftGPU.\n");
        throw "Cannot create SiftGPU.";
    }
}

unordered_map<size_t, pair<size_t, Vector2d>> SiftTracker::track(vector<unsigned char> &img) {
    unordered_map<size_t, pair<size_t, Vector2d>> result;

    static vector<SiftGPU::SiftKeypoint> oldKeys;
    static vector<array<float, 128>> oldDescriptors;

    m_sift->RunSIFT(m_img_width, m_img_height, img.data(), GL_LUMINANCE, GL_UNSIGNED_BYTE);
    int nfeat = m_sift->GetFeatureNum();
    vector<SiftGPU::SiftKeypoint> newKeys(nfeat);
    vector<array<float, 128>> newDescriptors(nfeat);

    m_sift->GetFeatureVector(newKeys.data(), newDescriptors[0].data());
    for (size_t i = 0; i < nfeat; ++i) {
        result[i].first = i;
        Vector2d p((newKeys[i].x - m_camera_cx) / m_camera_fx, (newKeys[i].y - m_camera_cy) / m_camera_fy);
        result[i].second = p;
    }

    if (newKeys.size() > 0 && oldKeys.size() > 0) {
        int matches[4096][2];

        m_matcher->SetDescriptors(0, (int)oldKeys.size(), oldDescriptors[0].data());
        m_matcher->SetDescriptors(1, (int)newKeys.size(), newDescriptors[0].data());

        int nmatch = m_matcher->GetSiftMatch(4096, matches);

        for (int i = 0; i < nmatch; ++i) {
            const SiftGPU::SiftKeypoint& oldKey = oldKeys[matches[i][0]];
            const SiftGPU::SiftKeypoint& newKey = newKeys[matches[i][1]];
            double dx = oldKey.x - newKey.x;
            double dy = oldKey.y - newKey.y;
            if (dx*dx + dy*dy < (m_img_width*m_img_height / 36)) {
                result[matches[i][1]].first = matches[i][0];
            }
        }

    }

    oldKeys.swap(newKeys);
    oldDescriptors.swap(newDescriptors);

    return result;
}
