% Comparing Adaptive Resolution and Fixed Resolution
%% select data 
datano = 4; % select data 1~4
Dtheta = -1.396394; R = 0.1803429337; % hardware parameters
switch(datano)
    case 1
        dir = 'data1\fixed\';
        dir_adp = 'data1\adaptive\';
        dataname = 'batchScanned20200703172752';
        dataname_adp = 'batchScanned20200703174420';
    case 2
        dir = 'data2\fixed\';
        dir_adp = 'data2\adaptive\';
        dataname = 'batchScanned20200703180037';
        dataname_adp = 'batchScanned20200703181723';
    case 3
        dir = 'data3\fixed\';
        dir_adp = 'data3\adaptive\';
        dataname = 'batchScanned20200703191706';
        dataname_adp = 'batchScanned20200703193139';
    case 4
        dir = 'data4\fixed\';
        dir_adp = 'data4\adaptive\';
        dataname = 'batchScanned20201024132319';
        dataname_adp = 'batchScanned20201024133958';
    otherwise % put new data here
        dir = 'data1\fixed\';
        dir_adp = 'data1\adaptive\';
        dataname = 'batchScanned20200703172752';
        dataname_adp = 'batchScanned20200703174420';
end

%% read data
[rData_raw, angleV, angleH_raw, timestamp] = ...
    read_scandata([dir,dataname,'.txt']); 
philen = length(angleH_raw);
[rData_raw_adp, angleV_adp, angleH_raw_adp, timestamp_adp] = ...
    read_scandata([dir_adp,dataname_adp,'.txt']); 
philen_adp = length(angleH_raw_adp);

Kdr = pi/180; % for converting degree to radian
Dtheta_s0 = 360/1024; % vertical angular resolution (deg)
Dphi_s0 = (angleH_raw(end)-angleH_raw(1))/philen;

%% show data
ps = range2points(rData_raw, angleV, angleH_raw, R, Dtheta);
figure(10);
scatter3(ps(:,1),ps(:,2),ps(:,3),1,'.'); 
xlim([-2,2]); ylim([-1.5 1]); zlim([-0.5,2])
az = 120; el = 80; view(az,el)
xlabel('x'); ylabel('y'); zlabel('z');
title([num2str(datano), ' fixed'])

% ps_adp = range2points(rData_raw_adp, angleV_adp, angleH_raw_adp, R, Dtheta);
% figure(11);
% scatter3(ps_adp(:,1),ps_adp(:,2),ps_adp(:,3),1,'.'); 
% xlim([-2,2]); ylim([-1.5 1]); zlim([-0.5,2])
% az = 120; el = 80; view(az,el)
% xlabel('x'); ylabel('y'); zlabel('z');
% title([num2str(datano), ' adaptive'])

dlen = philen_adp - philen;
figure(12); plot([[angleH_raw; zeros(dlen,1)] angleH_raw_adp]); 
title(['phi - ', num2str(datano)'])

%% calculate density
P = 600; % number of angle point
xL_max = 2; % max distance for searching objects
zL_max = 2; % max height for searching objects
dist_s = 1; % nominal working distance (m)
msgOff = false;

[pntSurf,surfPntSpc,surfDens,surfBeta,phi_p] = ...
    scanSurfaces(P,xL_max,zL_max,dist_s,msgOff, rData_raw,angleV,angleH_raw,R,Dtheta); 
[pntSurf_adp,surfPntSpc_adp,surfDens_adp,surfBeta_adp,phi_p_adp] = ...
    scanSurfaces(P,xL_max,zL_max,dist_s,msgOff, rData_raw_adp,angleV_adp,angleH_raw_adp,R,Dtheta); 

%% show results
% horizontal space between points
e_hori = dist_s*Dphi_s0*Kdr; % expected horizontal space between points
figure(1); 
plot(phi_p,surfPntSpc(:,1),...
    phi_p_adp,surfPntSpc_adp(:,1),...
    phi_p,e_hori*ones(P,1),'--k'); 
title('horizontal space between points (front)')
legend('fixed res.', 'adaptive res.', 'expected space')
figure(2); 
plot(phi_p,surfPntSpc(:,2),...
    phi_p_adp,surfPntSpc_adp(:,2),...
    phi_p,e_hori*ones(P,1),'--k'); 
title('horizontal space between points (back)')
legend('fixed res.', 'adaptive res.', 'expected space')

% density
figure(3); plot(phi_p, surfDens); 
title('point density (fixed res.)')
legend('positive plane', 'negative plane')
figure(4); plot(phi_p_adp, surfDens_adp); 
title('point density (adaptive res.)')
legend('positive plane', 'negative plane')

% surface position
figure(5);
scatter3(pntSurf(:,1),pntSurf(:,2),pntSurf(:,3),100,'.'); 
az = 22; el = 90; view(az,el)
xlabel('x'); ylabel('y'); zlabel('z');
title('detected object surfaces (fixed res.)')

figure(6);
scatter3(pntSurf_adp(:,1),pntSurf_adp(:,2),pntSurf_adp(:,3),100,'.'); 
az = 22; el = 90; view(az,el)
xlabel('x'); ylabel('y'); zlabel('z');
title('detected object surfaces (adaptive res.)')