- MSCKF

      包含了 A Multi-State Constraint Kalman Filter for Vision-aided Inertial Navigation 的核心算法实现

- MSCKF_Simulation

      使用模拟的数据来测试 MSCKF

- 坐标系

      设备坐标系（Gyro 和 Accelerometer）朝向
      将 iPhone 竖直放置，面对屏幕

           y
           |
           |
           +-----+
           |     |
           |     |
           |     |
           |  o  |
           +-----+---- x
          /
         /
        z

      相机坐标系朝向
      将 iPhone 横向放置，面对屏幕，Home键位于右侧

           z
          /
         /
        +--------------+---- x
        |              |
        |            O |
        |              |
        +--------------+
        |
        |
        y

      R_imu_to_img = [ -UnitY(), -UnitX(), -UnitZ() ];
      p_img_in_imu = [ 0.0065, 0.0638, 0.0000 ]; // 只是目测

- Accelerometer 含义

      “芯片对电路板施加的加速度（力/芯片质量）”
      因此当 iPhone 平放时，芯片会对电路板施加向下的力，对应加速度为设备坐标系的 -z 方向
      当 iPhone 向设备坐标系某个方向加速时，芯片对电路板施加反向的力，所以读数需要取负号