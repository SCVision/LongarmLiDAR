% Description: Let LIDAR system parallel to calibration board, using method of 
%              line fitting to calibrate theta.
% 
% Writen by LI, Shuqing (jay_issaac@163.com), 20191115
%

% clc;
% close all;
% clear;
path(path,'..\funcs')
path(path,'.\calibr_funcs')
%% data reading 
[range, angleV, angleH, timestamp] = read_scandata("scenario2_theta\batchScanned20191114103039.txt");
[n,m] = size(range);
%% display raw pointcloud 
R = 0.18; 
Dtheta = -2.6;
figure(20); show_pointcloud(range, angleV, angleH, R, Dtheta);
zlim([-0.1,3])
%% set plane boundary
range_denoise = remove_min_outlier(range, 0);
phi1 = 1; phi2 = 41;
theta1 = 292; theta2 = 428;
range_crop = range_denoise(phi1:phi2,theta1:theta2);
angleV_crop = angleV(theta1:theta2);
angleH_crop = angleH(phi1:phi2);

% figure(1); 
% plot(angleH); title('phi');
% figure(2); 
% plot(range'); title('range display');
% figure(3); 
% plot(range_crop'); title('range display');
% figure(21); show_pointcloud(range_crop, angleV_crop, angleH_crop, R, Dtheta);
% return; % comment the figures and the return after choose the scope
%% line fiting zL=k*xL+b 
d = mean(range_crop,1);
theta = angleV_crop;
[Dtheta,b]=calibrate_Dtheta(d, theta);
%% show fitting accuracy 
[nn,mm] = size(range_crop);
x = zeros(nn,mm);
z = zeros(nn,mm);
for i = 1:1:nn
   for j = 1:1:mm
       x(i,j) = range_crop(i,j)*sind(angleV_crop(j)+Dtheta);
       z(i,j) = range_crop(i,j)*cosd(angleV_crop(j)+Dtheta);       
   end
end
X = reshape(x.',nn*mm,1);
Z = reshape(z.',nn*mm,1);
Z1=zeros(size(X)) + b*cosd(Dtheta);
figure(4);
plot(X,[Z Z1],'.')
% ylim([0 3])
title('fitting accuracy');
error = Z1-Z;
mse = sqrt(mean(error.*error));
fprintf('MSE = %fm\n', mse)
%% display pointcloud 
% figure(22); show_pointcloud(range, angleV, angleH, R, Dtheta);