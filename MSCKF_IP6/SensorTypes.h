#pragma once

#include <cstdint>
#include <array>
#include <vector>
#include <Eigen/Eigen>
#include "UDPSocket.h"

#define SYNC_BUFFER_LENGTH 8
#define M_G 9.81

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

struct SensorSocketData {
    char event;
    double timestamp;
    Eigen::Vector3d p;
    double u;
    double v;
};
#pragma pack(pop)

class SyncBuffer {
public:
    void clear() {
        m_buffer.clear();
    }

    void push(const SensorData &data) {
        m_buffer.push_back(data);
        std::push_heap(m_buffer.begin(), m_buffer.end(), std::greater<SensorData>());
    }

    size_t size() const {
        return m_buffer.size();
    }

    bool canpop() const {
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

class DataFile {
public:
    DataFile(const std::string &filepath) {
        m_file = fopen(filepath.c_str(), "rb");
    }
    ~DataFile() {
        fclose(m_file);
    }

    bool eof() const {
        return feof(m_file) != 0;
    }

    bool get(SensorData &data) {
        if (eof()) return false;
        fread(&data, sizeof(SensorDataHeader), 1, m_file);
        if (eof()) return false;
        switch (data.event) {
        case EVENT_GPS:
            fread(&data.p, sizeof(Eigen::Vector3d), 1, m_file);
            fread(&data.u, sizeof(double), 1, m_file);
            fread(&data.v, sizeof(double), 1, m_file);
            break;
        case EVENT_ACCELEROMETER:
        {
            Eigen::Vector3d a;
            fread(&a, sizeof(Eigen::Vector3d), 1, m_file);
            data.p = -M_G*a;
        }
            break;
        case EVENT_GYRO:
            fread(&data.p, sizeof(Eigen::Vector3d), 1, m_file);
            break;
        case EVENT_MAGNETOMETER:
            fread(&data.p, sizeof(Eigen::Vector3d), 1, m_file);
            break;
        case EVENT_DEVICEMOTION:
            fread(&data.p, sizeof(Eigen::Vector3d), 1, m_file);
            break;
        case EVENT_CAMERA:
            data.img.resize(640 * 480);
            fread(data.img.data(), sizeof(unsigned char), 640 * 480, m_file);
            break;
        }
        return true;
    }
private:
    FILE *m_file = nullptr;
};

class DataSocket {
public:
    DataSocket(const std::string &notused) {
        m_socket.bind(udp::address(9123));
    }

    bool eof() const { return false; }

    bool get(SensorData &data) {
        udp::address from;
        while (m_socket.recv(from, &data, sizeof(SensorSocketData)) != sizeof(SensorSocketData));
        if (data.event == EVENT_ACCELEROMETER) {
            data.p = -M_G*data.p;
        }
        return true;
    }
private:
    udp::socket m_socket;
};

template<typename DataSource>
class DataReader {
public:
    DataReader(const std::string &filepath) : m_source(filepath) {}

    bool hasNext() {
        while (!m_sync.canpop()) {
            SensorData data;
            if (m_source.get(data)) {
                m_sync.push(data);
            }
            else {
                break;
            }
        }
        return m_sync.size() > 0;
    }
    SensorData next() {
        static double time = 0;
        SensorData data = m_sync.pop();
        if (data.timestamp < time) {
            std::cout << "Warning: timestamp interleaving!" << std::endl;
        }
        time = data.timestamp;
        return data;
    }

private:
    DataSource m_source;
    SyncBuffer m_sync;
};