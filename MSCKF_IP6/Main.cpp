#include <memory>
#include <Eigen/Eigen>
#include <opencv2/opencv.hpp>
#define WIN32_LEAN_AND_MEAN
#include <Windows.h>
#include <gl/GL.h>
#include <SiftGPU/SiftGPU.h>
#include <MSCKF.h>
#include "SensorTypes.h"

using namespace std;
using namespace Eigen;

bool nextData(FILE *file, SensorData &data) {
    if (feof(file)) {
        return false;
    }

    fread(&data, sizeof(SensorDataHeader), 1, file);
    if (data.event != EVENT_CAMERA) {
        fread(&data.p, sizeof(Vector3d), 1, file);
        if (data.event == EVENT_GPS) {
            fread(&data.u, sizeof(double), 2, file);
        }
    }
    else {
        vector<unsigned char> tmp;
        tmp.resize(640 * 480);
        data.img.resize(640 * 480);
        fread(tmp.data(), sizeof(unsigned char), 640 * 480, file);
        for (int y = 0; y < 480; ++y) {
            for (int x = 0; x < 640; ++x) {
                data.img[(x+1)*480-y-1] = tmp[y*640+x];
            }
        }
    }

    return true;
}

void imshow(const string &winname, const vector<unsigned char> &img) {
    cv::Mat i(640, 480, CV_8UC1, (void*)img.data());
    cv::imshow(winname, i);
}

double fx = 613.969482;
double fy = 613.969482;
double cx = 240.0;
double cy = 320.0;

void imshow(const string &winname, const vector<unsigned char> &img, const unordered_map<size_t, pair<size_t, Vector2d>> &match) {
    cv::Mat i(640, 480, CV_8UC1, (void*)img.data());
    for (size_t j = 0; j < match.size(); ++j) {
        const Vector2d &p = match.at(j).second;
        cv::circle(i, cv::Point(fx*p.x() + cx, fy*p.y() + cy), 3, cv::Scalar::all(0));
    }
    cv::imshow(winname, i);
}

class SiftTracker {
public:
    SiftTracker() {
        m_img_width = 480;
        m_img_height = 640;
        m_camera_fx = 613.969482;
        m_camera_fy = 613.969482;
        m_camera_cx = 240.0;
        m_camera_cy = 320.0;

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

    unordered_map<size_t, pair<size_t, Vector2d>> track(vector<unsigned char> &img) {
        static vector<SiftGPU::SiftKeypoint> oldKeys;
        static vector<array<float, 128>> oldDescriptors;

        unordered_map<size_t, pair<size_t, Vector2d>> result;

        m_sift->RunSIFT(m_img_width, m_img_height, img.data(), GL_LUMINANCE, GL_UNSIGNED_BYTE);
        int nfeat = m_sift->GetFeatureNum();
        vector<SiftGPU::SiftKeypoint> newKeys(nfeat);
        vector<array<float, 128>> newDescriptors(nfeat);

        m_sift->GetFeatureVector(newKeys.data(), newDescriptors[0].data());
        for (size_t i = 0; i < nfeat; ++i) {
            result[i].first = i;
            Vector2d p((newKeys[i].x-m_camera_cx)/m_camera_fx, (newKeys[i].y-m_camera_cy)/m_camera_fy);
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

private:
    int m_img_width;
    int m_img_height;

    double m_camera_fx;
    double m_camera_fy;
    double m_camera_cx;
    double m_camera_cy;

    unique_ptr<SiftGPU> m_sift;
    unique_ptr<SiftMatchGPU> m_matcher;
};

int main(int argc, char *argv[]) {
    FILE *file = nullptr;
    file = fopen(argv[1], "rb");

    SiftTracker tracker;

    SensorData data;
    while (nextData(file, data)) {
        if (data.event == EVENT_CAMERA) {
            auto matches = tracker.track(data.img);
            imshow("test", data.img, matches);
            cv::waitKey(1000 / 60);
        }
    }

    fclose(file);
    return 0;
}
