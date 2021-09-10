# LongarmLIDAR
Please cite this paper when using the data and codes:

[1] Jingyu Lin, Shuqing Li, Wen Dong, Takafumi Matsumaru, Shengli Xie. Long-arm three-dimensional LiDAR for anti-occlusion and anti-sparsity point clouds. IEEE Transactions on Instrumentation and Measurement, vol. 70, pp. 1-10, art no. 4506610, 2021. DOI: 10.1109/TIM.2021.3104019

Data and codes for long-arm 3D LIDAR

- Demos: Read and show LIDAR data. 

- Calibration: Process of calibration with two degree of freedom. Run 'main_calibratetheta2' first and then 'main_calibrateR'.

- CalibrationPrecise: Process of calibration with four degree of freedom. Run 'main_calibrateAll_201357' first and then 'main_calibrate_byhand201357'.

- Anti_occlusion_test: Data of anti-occlusion scenarios and display functions. Objects are box, sphere and column, distances are from 350mm to 4000mm. 

- Adaptive_scanning: Demonstration for adaptive resolution scanning. Objects are box, sphere and column, distances are from 600mm to 2000mm. 

- Funcs: General functions supporting commands in the other paths.

- LiDAR_sim: Simulate adaptive resulition scanning or fixed resulition scanning of a long-arm lidar in a 2D map (top view).
