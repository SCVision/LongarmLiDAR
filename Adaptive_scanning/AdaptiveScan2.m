% Adaptive Resolution Scanning for Long-arm LIDAR
% This is a demonstration program. Input a scanned data (usually scanned 
% with a fixed step), and it will do the following work:
% 1) update angular steps Dphi_s based on point density
% 2) process Dphi_s and generate executable scanning commands for LIDAR
% 3) show new angular steps and the positions of objects.
%
% Key parameters:
% R, Dtheta - calibrated hardware parameters
% dataname - data file name
% nPhi1, nPhi2 - used to trim head and tail of data
% dist_s - working distance
% xL_max, zL_max - max distance and height for searching objects
% P - numbers for checking density and updating angular steps
% W_smooth, W_seg - control the effect of step processing
%
% Output results
% Dphi_s - orginal new angular steps
% Dphi_s_final - processed new angular steps
% pntSurf - point set of the positions of detected surfaces
%
% Writen by LIN, Jingyu (linjy02@hotmail.com), 20200622
%
path(path,'..\Funcs')
path(path,'.\adaptive_funcs')

%% select data 
datano = 1; % select data 1~3
bReverse = false; % virtually scan from the other end
Dtheta = -1.396394; R = 0.1803429337; % hardware parameters
switch(datano)
    case 1
        dir = 'data1\fixed\';
        dataname = 'batchScanned20200703172752';
        nPhi1 = 21; nPhi2 = 3160; % choose horizontal data range
    case 2
        dir = 'data2\fixed\';
        dataname = 'batchScanned20200703180037';
        nPhi1 = 21; nPhi2 = 3160; % choose horizontal data range
    case 3
        dir = 'data3\fixed\';
        dataname = 'batchScanned20200703191706';
        nPhi1 = 21; nPhi2 = 3160; % choose horizontal data range
    otherwise % put new data here
        dir = 'data1\fixed\';
        dataname = 'batchScanned20200703172752';
        nPhi1 = 21; nPhi2 = 3160; % choose horizontal data range
end

% read data
[range_raw, angleV, angleH_raw, timestamp] = ...
    read_scandata([dir,dataname,'.txt']); 

% whether virtually scan from the other end  
if bReverse
    range_raw = range_raw(end:-1:1, :);
    angleH_raw = angleH_raw(end:-1:1);
end
range = range_raw;
angleH = angleH_raw;

% show data
ps = range2points(range_raw, angleV, angleH_raw, R, Dtheta);
figure(10);
scatter3(ps(:,1),ps(:,2),ps(:,3),1,'.'); 
xlim([-2,2]); ylim([-1.5 1]); zlim([-0.5,2])
az = 120; el = 80; view(az,el)
xlabel('x'); ylabel('y'); zlabel('z');
title(dataname)
% figure(11); plot(angleH_raw); title(['phi - ', dataname])
% return

%% setup data
NPhi = nPhi2 - nPhi1 + 1; % number of horizontal data points
phi_step = (angleH(nPhi2)-angleH(nPhi1))/NPhi; % horizontal angular step (deg)

% make sure angleH start from 0 and move clockwise
bCounterClock = (phi_step<0); % whether scan counter clockwise
if bCounterClock
    % setup a mirror environment
    phi_step = -phi_step;
    angleH = -angleH;
end
angleH0 = angleH(1);
angleH = angleH - angleH0;

%% setup data and parameters
% 2. constants
Kdr = pi/180; % for converting degree to radian
Dtheta_s0 = 360/1024; % vertical angular resolution (deg)

% 3. parameters
dist_s = 1; % nominal working distance (m)
Dphi_s0 = phi_step; % horizontal angular step (deg)
e_hori = dist_s*Dphi_s0*Kdr; % expected horizontal space between points
dens_s = 1/(e_hori*dist_s*Dtheta_s0*Kdr); % expected density
hG = 6*e_hori; % width of Gaussian (m)

xL_max = 2; % max distance for searching objects
zL_max = 2; % max height for searching objects
xL = 0.4:0.02:xL_max; % horizontal points in scanning plane for density
zL = 0:0.1:zL_max; % vertical points in scanning plane for density

Dphi_rho = 4; % phi range for calculating density (deg)
% Dphi_rho = 30*phi_step; % phi range for calculating density (deg)
NDphi_rho = ceil(Dphi_rho / phi_step); % number of data for calculating density

P = 600; % number of rotation positions

%% initialize arrays
phi_p_step = 360/P;
phi_p = (1:P)'*phi_p_step; % rotation positions
Dphi_s = Dphi_s0*ones(P,1); % record horizontal angular step (deg)
surfPos_p = zeros(P,2); % record surface at + scanning plane (xL, yL)
surfPos_n = zeros(P,2); % record surface at - scanning plane (xL, yL)
surfDens = zeros(P,2); % record point density on surface at scanning plane (+ & -)

%% calculate new angular step
Ntheta = length(angleV); % number of vertical data points
rangV_p = Ntheta/2:Ntheta; % vertical angle for position direction
rangV_n = 1:Ntheta/2;      % vertical angle for negative direction

fprintf('Updating angular steps. \n');
n_step = ceil(phi_p_step/phi_step); 
for n = nPhi1 : n_step: nPhi2
    rangH = n+(0:NDphi_rho);
    phi = (angleH(n) + angleH(n+NDphi_rho))/2; % orientation
    idx_p = 1 + round(phi / phi_p_step);
    
    % positive orientation, +xL, +angleV
    x_p = xL*cosd(phi) + R*sind(phi);
    y_p = -xL*sind(phi) + R*cosd(phi);
    ps_crop = range2points(range(rangH,rangV_p), angleV(rangV_p), angleH(rangH), R, Dtheta);
%     figure(11);
%     scatter3(ps_crop(:,1),ps_crop(:,2),ps_crop(:,3),1,'.');
%     xlabel('x'); ylabel('y'); zlabel('z');
    pln_dens = KDEonVplane([x_p;y_p]', zL, ps_crop, hG); 

    % update angular resolution Dphi_s
    [dm, zm, objDens, beta] = findSurface(pln_dens, xL, zL, Dphi_s0, Dtheta_s0); 
    surfDens(idx_p,1) = objDens;
    surfPos_p(idx_p,1) = dm;
    surfPos_p(idx_p,2) = zm;
    if dm  > R
        alpha = atand(dm/R);
        phi_q = phi + 2*alpha; % positive direction
        idx_q = round(phi_q / phi_p_step);
%         if phi_q < 360
        if idx_q <= P
            Dphi_sq = abs(e_hori * sind(2*alpha-beta) / dm) / Kdr; % degree
            Dphi_s(idx_q) = min(Dphi_s(idx_q), Dphi_sq);
            
            % output message
            fprintf('p=%.1f d%.0fk +q %.1f, beta %.1f -- %.1f, step %.4f\n',...
                phi,objDens/1000,phi_q,beta,beta+180-2*alpha,Dphi_sq );
        else
            % output message
            fprintf('p=%.1f d%.0fk +q out of range\n', phi,objDens/1000);
        end
    else
        % output message
        fprintf('p=%.1f d%.0fk +q no object\n', phi,objDens/1000);
    end
    
    % negative orientation, -xL, -angleV
    x_n = -xL*cosd(phi) + R*sind(phi);
    y_n = xL*sind(phi) + R*cosd(phi);
    ps_crop = range2points(range(rangH,rangV_n), angleV(rangV_n), angleH(rangH), R, Dtheta);
%     figure(12);
%     scatter3(ps_crop(:,1),ps_crop(:,2),ps_crop(:,3),1,'.');
%     xlabel('x'); ylabel('y'); zlabel('z');
    pln_dens = KDEonVplane([x_n;y_n]', zL, ps_crop, hG); 

    % update angular resolution Dphi_s
    [dm, zm, objDens, beta] = findSurface(pln_dens, xL, zL, Dphi_s0, Dtheta_s0); 
    surfDens(idx_p,2) = objDens;
    surfPos_n(idx_p,1) = dm;
    surfPos_n(idx_p,2) = zm;
    if dm  > R
        alpha = atand(dm/R);
        phi_q = phi - 2*alpha; % negative direction 
        phi_q = phi_q + 360; % phi_q must be negative
        idx_q = round(phi_q / phi_p_step);
%         if phi_q<0
        if idx_q <= P
            Dphi_sq = abs(e_hori * sind(2*alpha-beta) / dm) / Kdr; % degree
            Dphi_s(idx_q) = min(Dphi_s(idx_q), Dphi_sq);
            
            % output message
            fprintf('p=%.1f d%.0fk -q %.1f, beta %.1f -- %.1f, step %.4f\n',...
                phi,objDens/1000,phi_q,beta,beta+180-2*alpha,Dphi_sq );
        else
            % output message
            fprintf('p=%.1f d%.0fk -q out of range\n', phi,objDens/1000);
        end
    else
        % output message
        fprintf('p=%.1f d%.0fk -q no object\n', phi,objDens/1000);
    end
end
fprintf('End of updating. \n\n');

%% process new angular steps
W_smooth = 3; % control smooth effect
n = (-W_smooth:W_smooth);
Dphi_s_smooth = Dphi_s;
for i=1+W_smooth:length(Dphi_s)-W_smooth
    m = Dphi_s(i+n);
    Dphi_s_smooth(i) = min(m);
end

% segment
W_seg = 16; % control segmentation effect
thr = 0.05; % control segmentation
n = (1:W_seg);
Dphi_s_final = Dphi_s_smooth;
D_Dphi_s = abs(diff(Dphi_s_smooth));
diffmax = max(D_Dphi_s);
LIDAR_T0 = 0.1; % interval of LIDAR scanning
np0 = 1;
fprintf('Scanning commands: \nstart, end, velocity\n');
for i=1+W_seg : W_seg : length(Dphi_s)-W_seg
    [m, idx] = max(D_Dphi_s(i+n));
    np1 = i+idx;
    if m > thr*diffmax
        % find a segment
        avgDphiSeg = mean(Dphi_s_smooth(np0:np1));
        Dphi_s_final(np0:np1) = avgDphiSeg;
        fprintf('%.2f, %.2f, %.2f\n',phi_p(np0), phi_p(np1+1), avgDphiSeg/LIDAR_T0);
        np0 = np1 + 1;
    end
end
np1 = length(Dphi_s);
avgDphiSeg = mean(Dphi_s_smooth(np0:np1));
fprintf('%.1f, %.1f, %.1f\n',phi_p(np0), phi_p(np1), avgDphiSeg/LIDAR_T0);
Dphi_s_final(np0:end) = avgDphiSeg;
fprintf('End of scanning commands\n');

figure(1); plot(phi_p, [Dphi_s Dphi_s_smooth],phi_p,Dphi_s_final,'k'); 
ylim([0, Dphi_s0*1.1]);
title('adaptive angular steps'); 
legend('origin','smoothed','final')

%% show results
% density
figure(2); plot(phi_p, [surfDens dens_s*ones(P,1)]); 
title('point density at object surface')
legend('positive plane', 'negative plane', 'expected density')

% surface position
pntSurf_p = zeros(P,3); % positions at + scanning plane
pntSurf_n = zeros(P,3); % positions at - scanning plane
SH = sind(phi_p); 
CH = cosd(phi_p); 
for i = 1:P
    if surfPos_p(i,1) > 0.02
        pntSurf_p(i,1) =  surfPos_p(i,1)*CH(i)+ R*SH(i);
        pntSurf_p(i,2) = -surfPos_p(i,1)*SH(i)+ R*CH(i);
        pntSurf_p(i,3) =  surfPos_p(i,2);
    end
    if surfPos_n(i,1) > 0.02
        pntSurf_n(i,1) = -surfPos_n(i,1)*CH(i)+ R*SH(i);
        pntSurf_n(i,2) =  surfPos_n(i,1)*SH(i)+ R*CH(i);
        pntSurf_n(i,3) =  surfPos_n(i,2);
    end
end
figure(3);
pntSurf = [ pntSurf_p; pntSurf_n];
scatter3(pntSurf(:,1),pntSurf(:,2),pntSurf(:,3),100,'.'); 
az = 22; el = 90; view(az,el)
xlabel('x'); ylabel('y'); zlabel('z');
title('detected object surfaces')
