function [pnts_p,pnts_n] = range2points_topview(pos, phi, ...
    range_p, range_n, arm_len) 
% Function: convert long-arm LIDAR range data (topview) to 2D point sets.
% Input:
%     pos - X, Y coordinates of the LiDAR in the map.
%     phi - array of long-arm orientaion (Nx1), deg.
%     range_p, range_n - array of range corresponding to phi (Nx1).
%     arm_len - arm length.
% Output:
%     pnts_p,pnts_n - X, Y coordinates of points (N*2)
% Demo:
% [rData, angleV, angleH, timestamp] = read_scandata('Scanned1.txt'); 
% iPhi = 100;
% pnts = range2points2D(rData(iPhi,:)', angleV);
% figure(13); plot(pnts(:,1), pnts(:,2),'b.')
% 
% Writen by LIN, Jingyu (linjy02@hotmail.com), 20210119
%
SV = sind(phi(:)); CV = cosd(phi(:)); 
xL = pos(1) + arm_len.*CV; % pos of 2D LiDAR
yL = pos(2) - arm_len.*SV; 
x = xL + range_p.*-SV; % positive direction: phi+90
y = yL - range_p.*CV;
pnts_p = [x(:), y(:)];
x = xL + range_n.*SV; % negative direction: phi-90
y = yL - range_n.*-CV;
pnts_n = [x(:), y(:)];
