% Description: calibration of R, Dphi, Dpsi, Dtheta
% Method: use z-x-y nautical angles 
% Steps: 
% 1. find a ceiling and a wall.
% 2. scan from 0 deg to 360 deg.
% 3. set 'cropping=true', crop calibration planes on wall and ceiling.
% 4. set 'cropping=false', run the code.
% 
path(path,'..\funcs')
cropping = false;

%% data reading 
[range, angleV, angleH, timestamp] = ...
    read_scandata("calibr_data\batchScanned20200716201357.txt");
%% display raw pointcloud 
R = 0.2;
Dphi = 0;
Dpsi = 0;
Dtheta = 0;
ps0 = range2pointsPrecise(range, angleV, angleH, R,Dphi,Dpsi,Dtheta);
figure(1); scatter3(ps0(:,1),ps0(:,2),ps0(:,3),1,'.'); 
xlabel('x'); ylabel('y'); zlabel('z');
% figure(2); plot(angleH); 
minRange = 0.1;

%% crop calibration plane
% wall
phi1_1 = 1540; phi1_2 = 1850; theta1_1 = 600;  theta1_2 = 640;
phi2_1 = 3070; phi2_2 = 3200; theta2_1 = 53; theta2_2 = 90;
phi3_1 = 21; phi3_2 = 140; theta3_1 = 53; theta3_2 = 90;
range_crop1 = range(phi1_1:phi1_2,theta1_1:theta1_2);
angleH_crop1 = angleH(phi1_1:phi1_2);
angleV_crop1 = angleV(theta1_1:theta1_2);
range_crop2 = range(phi2_1:phi2_2,theta2_1:theta2_2);
angleH_crop2 = angleH(phi2_1:phi2_2);
angleV_crop2 = angleV(theta2_1:theta2_2);
range_crop3 = range(phi3_1:phi3_2,theta3_1:theta3_2);
angleH_crop3 = angleH(phi3_1:phi3_2);
angleV_crop3 = angleV(theta3_1:theta3_2);

% ceiling
% phi5_1 = 1700; phi5_2 = 2900; theta5_1 = 318; theta5_2 = 400;
phi6_1 = 21; phi6_2 = 3220; theta6_1 = 290; theta6_2 = 410;
c_theta = length(angleV)/2; c_remove = 20;
theta_domain = [theta6_1:(c_theta-c_remove),(c_theta+c_remove):theta6_2];
range_crop6 = range(phi6_1:phi6_2,theta_domain);
angleH_crop6 = angleH(phi6_1:phi6_2);
angleV_crop6 = angleV(theta_domain);

% % remove outlier
% range_crop6 = remove_outlier(range_crop6,minRange);

% show patches
if cropping
    % show wall patches
    ps1 = range2pointsPrecise(range_crop1, angleV_crop1, angleH_crop1, R,Dphi,Dpsi,Dtheta);
    ps2 = range2pointsPrecise(range_crop2, angleV_crop2, angleH_crop2, R,Dphi,Dpsi,Dtheta);
    ps3 = range2pointsPrecise(range_crop3, angleV_crop3, angleH_crop3, R,Dphi,Dpsi,Dtheta);
    ps_wall = [ps1;ps2;ps3];
    figure(3); scatter3(ps_wall(:,1),ps_wall(:,2),ps_wall(:,3),1,'.');
    xlabel('x'); ylabel('y'); zlabel('z');
    az = 90; el = 0; view(az,el)
    
    % show ceiling patches
    % ps5 = range2pointsPrecise(range_crop5, angleV_crop5, angleH_crop5, R,Dphi,Dpsi,Dtheta);
    ps6 = range2pointsPrecise(range_crop6, angleV_crop6, angleH_crop6, R,Dphi,Dpsi,Dtheta);
    ps_ceil = [ps6];
    figure(4); scatter3(ps_ceil(:,1),ps_ceil(:,2),ps_ceil(:,3),1,'.');
    xlabel('x'); ylabel('y'); zlabel('z');
    az = 0; el = 90; view(az,el)
    return; % comment the figures and the return after choose the scope
end

%% calibration Dtheta
% ceiling data
d_ceil = range_crop6;
theta_ceil = angleV_crop6;
phi_ceil = angleH_crop6;

% calibration
msgOn = true;
iter = 40;
Dtheta_step = 1;
[Dtheta, err1] = calibrate_Dtheta(d_ceil,theta_ceil,phi_ceil, ...
    R, Dphi, Dpsi, Dtheta, iter, Dtheta_step, msgOn) ;

%% calibration R & Dphi & Dpsi
% front data
d1 = range_crop1;
theta1 = angleV_crop1;
phi1 = angleH_crop1;

% back data
d2 = [range_crop2; range_crop3];
theta2 = angleV_crop2;
phi2 = [angleH_crop2; angleH_crop3];

iter = 60;
lr = 0.1;
[R, Dphi, Dpsi, err2] = calibrate_R_DphiDpsi(d1,theta1,phi1,d2,theta2,phi2, ...
    R, Dphi, Dpsi, Dtheta, iter, lr, msgOn) ;

%% display raw pointcloud 
ps0_cabr = range2pointsPrecise(range, angleV, angleH, R,Dphi,Dpsi,Dtheta);
figure(11); scatter3(ps0_cabr(:,1),ps0_cabr(:,2),ps0_cabr(:,3),1,'.'); 
xlabel('x'); ylabel('y'); zlabel('z');
% xlim([-2,2]); ylim([-1.5 1]); zlim([-0.5,2])
% az = 120; el = 80; view(az,el)

% show the patches
ps1_cabr = range2pointsPrecise(range_crop1, angleV_crop1, angleH_crop1, R,Dphi,Dpsi,Dtheta);
ps2_cabr = range2pointsPrecise(range_crop2, angleV_crop2, angleH_crop2, R,Dphi,Dpsi,Dtheta);
ps3_cabr = range2pointsPrecise(range_crop3, angleV_crop3, angleH_crop3, R,Dphi,Dpsi,Dtheta);
ps_cabr = [ps1_cabr;ps2_cabr;ps3_cabr];
figure(13); scatter3(ps_cabr(:,1),ps_cabr(:,2),ps_cabr(:,3),1,'.'); 
xlabel('x'); ylabel('y'); zlabel('z'); 
az = 90; el = 0; view(az,el)

ps_ceil_cabr = range2pointsPrecise(range_crop6, angleV_crop6, angleH_crop6, R,Dphi,Dpsi,Dtheta);
figure(14); scatter3(ps_ceil_cabr(:,1),ps_ceil_cabr(:,2),ps_ceil_cabr(:,3),1,'.'); 
xlabel('x'); ylabel('y'); zlabel('z'); 
az = 0; el = 90; view(az,el)
