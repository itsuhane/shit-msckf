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

- (11) 中 P_{IC} 和 P_{CC} 矩阵的结构

      设 J_c 为 (16) 中 J 矩阵左侧的 15 列，另外用 (V;k) 表示第 k 时间 V 的值
      那么根据 (15) 和 (16) 可以知道：

      P_{IC}|k 可以每 6 列分块，其中每块内容为 (P_{II};i)*(J_c;i)^T。
      P_{CC} 矩阵是一个对称矩阵，对它按 6x6 为单位分块的话，那么第 i 行 j 列的块为 (J_c;i)*(P_{II};min(i,j))*(J_c;j)^T。

      因此我们可以在第 i 个相机状态中伴随保存 (P_{II};i) 和 (J_c;i) 。在需要用到 P_{IC} 和 P_{CC} 的时候重构两个矩阵。

      （由于P_II和P_IC会随时间更新，这个方法不能直接使用，暂时搁置）