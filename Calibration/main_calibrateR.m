% Description: Calibrate R. 
% Setup: The LIDAR is horizontal and a calibration board is vertical.
%    Scan by left side and right side of the 2D LIDAR on the same plane
%    of Ax+By+Cz=1, respectively. 
% Writen by LIN Jingyu (linjy02@hotmail.com), 20200129
%

% clc;
% close all;
% clear;
path(path,'..\funcs')
path(path,'.\calibr_funcs')
%% data reading 
[range, angleV, angleH, timestamp] = read_scandata("scenario4_R\batchScanned20191114104048.txt");
%% display raw pointcloud 
% Dtheta = -2.677827; % accurate value
Dtheta = 0; % pay attention to the fraction due to inaccurate Dtheta
R = 0.2;
figure(20); show_pointcloud(range, angleV, angleH, R, Dtheta);
%% set plane boundary
range_denoise = remove_min_outlier(range, 0);
% range_denoise = range;
% phi1 = 1300; phi2 = 1540; theta1 = 130; theta2 = 188;
% phi1 = 20; phi2 = 340; theta1 = 510; theta2 = 560;
phi1_1 = 10; phi1_2 = 340; theta1_1 = 510; theta1_2 = 560;
range_crop1 = range_denoise(phi1_1:phi1_2,theta1_1:theta1_2);
angleH_crop1 = angleH(phi1_1:phi1_2);
angleV_crop1 = angleV(theta1_1:theta1_2);
phi2_1 = 1220; phi2_2 = 1550; theta2_1 =130;  theta2_2 = 190;
range_crop2 = range_denoise(phi2_1:phi2_2,theta2_1:theta2_2);
angleH_crop2 = angleH(phi2_1:phi2_2);
angleV_crop2 = angleV(theta2_1:theta2_2);

% figure(21); show_pointcloud(range_crop1, angleV_crop1, angleH_crop1, R, Dtheta);
% xlim([0,3]); zlim([0,3])
% figure(22); show_pointcloud(range_crop2, angleV_crop2, angleH_crop2, R, Dtheta);
% xlim([0,3]); zlim([0,3])
% return; % comment the figures and the return after choose the scope
%% plane fiting Ax+By+Cz=1 
iter = 1000;
d1 = range_crop1;
theta1 = angleV_crop1;
phi1 = angleH_crop1;
d2 = range_crop2;
theta2 = angleV_crop2;
phi2 = angleH_crop2;
[R,error] = calibrate_R(d1, theta1, phi1, d2, theta2, phi2, Dtheta, iter); 
fprintf('R = %.10f, error = %.10f\n', R, error(end));
%% display loss curve 
figure(1);
plot(error);
title('log loss curve');
%% display raw pointcloud 
Dtheta = -2.677827; % accurate value
R = 0.1750963106; % accurate value
figure(30); show_pointcloud(range, angleV, angleH, R, Dtheta);
% range_processed=medfilt2(range,[5 5]); % outlier removed
% figure(30); show_pointcloud(range_processed, angleV, angleH, R, Dtheta);

