% Description: use any data to manually calibrate R, Dphi, Dpsi, Dtheta 
% 1. select LIDAR data
% 2. setup range of ceiling and walls
% 3. adjust Dtheta, Dpsi to make the ceiling overlapped
% 4. adjust R, Dphi to make the walls overlapped
% Then repeat 3 & 4 until all objects perfectly overlapped
path(path,'..\Funcs')

%% 1. read data
dir = 'calibr_data';
dataname = 'batchScanned20200626210606.txt';
[range, angleV, angleH, timestamp] = read_scandata([dir,'\',dataname]); 

%% 2. initiallize parameters
ceil_z1 = 2; ceil_z2 = 3.5; % height of ceiling
wall_x1 = -3.5; wall_x2 = 4; % range of walls
wall_y1 = -3.5; wall_y2 = 4; % range of walls 
R=0.2; Dphi=0; Dpsi=0; Dtheta=0;

%% 3. manually calibrate Dtheta, Dpsi 
Dpsi=0.0269; Dtheta=-2.5;
ps = range2pointsPrecise(range, angleV, angleH, R,Dphi,Dpsi,Dtheta);
figure(20); scatter3(ps(:,1),ps(:,2),ps(:,3),1,'.'); 
az = 177; el = 2; view(az,el)
xlim([wall_x1,wall_x2]); ylim([wall_y1 wall_y2]); 
zlim([ceil_z1,ceil_z2])
% ylim([-2 2]); zlim([-.5,1])
xlabel('x'); ylabel('y'); zlabel('z');

%% 4. manually calibrate R, Dphi 
R=0.193821; Dphi=0.1598;
ps = range2pointsPrecise(range, angleV, angleH, R,Dphi,Dpsi,Dtheta);
figure(21); scatter3(ps(:,1),ps(:,2),ps(:,3),1,'.'); 
az = 180; el = 90; view(az,el)
xlabel('x'); ylabel('y'); zlabel('z');
xlim([wall_x1,wall_x2]); ylim([wall_y1 wall_y2]); 
zlim([-1, ceil_z1])
% zlim([-1, 1])

%% output parameters 
fprintf('R=%.6f; Dphi=%.4f; Dpsi=%.4f; Dtheta=%.4f;\n',R,Dphi, Dpsi, Dtheta)
