#include <opencv2/opencv.hpp>
#include <MSCKF.h>

#include "SiftTracker.h"
#include "AccelerometerCalibrator.h"

#include "../MSCKF_Simulation/Helpers.h"

#include "SensorTypes.h"

using namespace std;
using namespace Eigen;

int main(int argc, char *argv[]) {
    udp::socket::startup();
    DataReader<DataSocket> reader(argv[1]);

    while (reader.hasNext()) {
        SensorData data = reader.next();
    }

    udp::socket::cleanup();
    return 0;
}
