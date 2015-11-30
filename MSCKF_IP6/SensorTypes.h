#pragma once

#include <cstdint>
#include <array>
#include <vector>
#include <Eigen/Eigen>

#define SYNC_BUFFER_LENGTH 4

enum SensorEvent {
    EVENT_GPS = 1,
    EVENT_ACCELEROMETER = 2,
    EVENT_GYRO = 3,
    EVENT_MAGNETOMETER = 4,
    EVENT_DEVICEMOTION = 5,
    EVENT_CAMERA = 6
};

#pragma pack(push)
#pragma pack(1)
struct SensorDataHeader {
    char event;
    double timestamp;
};

struct RawSensorData {
    char event;
    double timestamp;
    Eigen::Vector3d p;
    double u;
    double v;
};

struct SensorData {
    char event;
    double timestamp;
    Eigen::Vector3d p;
    double u;
    double v;
    std::vector<unsigned char> img;

    bool operator< (const SensorData &d) const {
        if (timestamp != d.timestamp) {
            return timestamp < d.timestamp;
        }
        else {
            return event < d.event;
        }
    }
    bool operator>(const SensorData &d) const {
        return d < *this;
    }
};

struct ServerData {
    std::uint32_t id;
    char type;
    std::array<unsigned char, 3> color;
    Eigen::Vector4d vec0;
    Eigen::Vector4d vec1;
    Eigen::Vector4d vec2;
};
#pragma pack(pop)

class sync_buffer {
public:
    void clear() {
        m_buffer.clear();
    }

    void push(const SensorData &data) {
        m_buffer.push_back(data);
        std::push_heap(m_buffer.begin(), m_buffer.end(), std::greater<SensorData>());
    }

    bool can_pop() const {
        return m_buffer.size() > SYNC_BUFFER_LENGTH;
    }

    SensorData pop() {
        std::pop_heap(m_buffer.begin(), m_buffer.end(), std::greater<SensorData>());
        SensorData result = m_buffer.back();
        m_buffer.pop_back();
        return result;
    }

private:
    std::vector<SensorData> m_buffer;
};
