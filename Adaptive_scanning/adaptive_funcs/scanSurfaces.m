function [pntSurf,surfPntSpc,surfDens,surfBeta,phi_p] = scanSurfaces(P,xL_max,zL_max,dist_s,msgOff, rData,angleV,angleH,R,Dtheta) 
% Function: scan around and find object surfaces.
% Input:
%     P - number of rotation positions
%     xL_max, zL_max - max distance & height for searching objects
%     dist_s - nominal working distance
%     msgOff - true to close echo message
%     range - range data (H*V). 
%     angleV - vertical angles theta (V*1).
%     angleH - horizontal angles phi (H*1). 
%     R, Dtheta - calibration parameters 
% Output:
%     pntSurf - point set indicating position of surfaces (2P*3)
%     surfPntSpc - Hspace between points in +&- direction (P*2)
%     surfDens - points density in +&- direction (P*2)
%     surfBeta - surface orientation in +&- direction (P*2)
%     phi_p - rotation positions
% Test Demo: none
% 
% Writen by LIN, Jingyu (linjy02@hotmail.com), 20200707
%
Kdr = pi/180; % for converting degree to radian
Dtheta_s0 = 360/1024; % vertical angular resolution (deg)
Dphi_rho = 3; % phi range for calculating density (deg)

philen = length(angleH);
bCounterClock = (angleH(end) < angleH(1)); % whether scan counter clockwise
if bCounterClock
    % setup a mirror environment
    angleH = -angleH;
end
angleH = angleH - angleH(1); % phi starts from 0deg and is ascendent 

%% show data
% ps = range2points(range, angleV, angleH, R, Dtheta);
% figure(100);
% scatter3(ps(:,1),ps(:,2),ps(:,3),1,'.'); 
% xlim([-2,2]); ylim([-1.5 1]); zlim([-0.5,2])
% az = 120; el = 80; view(az,el)
% xlabel('x'); ylabel('y'); zlabel('z');

%% calculate density
% scanning plane
xL = 0.4:0.02:xL_max; % horizontal points in scanning plane for density
zL = 0:0.1:zL_max; % vertical points in scanning plane for density

% angle points
phi_p_step = (angleH(end))/P;
phi_p = (1:P)'*phi_p_step; % rotation positions
surfPos_p = zeros(P,2); % record surface at + scanning plane (xL, yL)
surfPos_n = zeros(P,2); % record surface at - scanning plane (xL, yL)
surfPntSpc = zeros(P,2); % record Hspace between points on surface at scanning plane (+ & -)
surfDens = zeros(P,2); % record point density on surface at scanning plane (+ & -)
surfBeta = zeros(P,2); % record orientation angle of surface at scanning plane (+ & -)

% calculate density
Ntheta = length(angleV); % number of vertical data points
rangV_p = Ntheta/2:Ntheta; % vertical angle for position direction
rangV_n = 1:Ntheta/2;      % vertical angle for negative direction
hG = 6*dist_s*(angleH(end)/philen)*Kdr; % width of Gaussian (m)
dataPC1 = 0;
if ~msgOff
    fprintf('Start searching surfaces.\n');
end
for n = 1 : P
    phi = phi_p(n);
    x_p = xL*cosd(phi) + R*sind(phi);
    y_p = -xL*sind(phi) + R*cosd(phi);
    x_n = -xL*cosd(phi) + R*sind(phi);
    y_n = xL*sind(phi) + R*cosd(phi);
    
    % get angle range for density
    for i = dataPC1+1 : philen
        if angleH(i) > phi-Dphi_rho
            break; 
        end
    end
    dataPC1 = i; 
    for i = dataPC1+1 : philen
        if angleH(i) > phi+Dphi_rho
            break; 
        end
    end
    dataPC2 = i-1; 
    rangH = dataPC1:dataPC2;
    Dphi_s0 = (angleH(dataPC2)-angleH(dataPC1))/length(rangH);
    
    % positive orientation, +xL, +angleV
    ps_crop = range2points(rData(rangH,rangV_p), ...
        angleV(rangV_p), angleH(rangH), R, Dtheta);
    pln_dens = KDEonVplane([x_p;y_p]', zL, ps_crop, hG); 

    % find object surface
    [dm, zm, objDens, beta] = findSurface(pln_dens, xL, zL, Dphi_s0, Dtheta_s0); 
    surfDens(n,1) = objDens;
    surfBeta(n,1) = beta;
    surfPos_p(n,1) = dm;
    surfPos_p(n,2) = zm;
    if dm>0 
        surfPntSpc(n,1) = 1/(objDens*dm*Dtheta_s0*Kdr);
        if ~msgOff
            fprintf('phi%5.1f +: ps %.1fmm, dens %2.0fk, (%.2f, %.2f)\n',...
                phi, surfPntSpc(n,1)*1000, objDens/1000, dm, zm);
        end
    end
    
    % negative orientation, -xL, -angleV
    ps_crop = range2points(rData(rangH,rangV_n),...
        angleV(rangV_n), angleH(rangH), R, Dtheta);
    pln_dens = KDEonVplane([x_n;y_n]', zL, ps_crop, hG); 

    % find object surface
    [dm, zm, objDens, beta] = findSurface(pln_dens, xL, zL, Dphi_s0, Dtheta_s0); 
    surfDens(n,2) = objDens;
    surfBeta(n,2) = beta;
    surfPos_n(n,1) = dm;
    surfPos_n(n,2) = zm;
    if dm>0
        surfPntSpc(n,2) = 1/(objDens*dm*Dtheta_s0*Kdr);
        if ~msgOff
            fprintf('phi%5.1f -: ps %.1fmm, dens %2.0fk, (%.2f, %.2f)\n',...
                phi, surfPntSpc(n,2)*1000, objDens/1000, dm, zm);
        end
    end
end
if ~msgOff
    fprintf('End of searching.\n');
end

%% transform surface position
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
pntSurf = [ pntSurf_p; pntSurf_n];
