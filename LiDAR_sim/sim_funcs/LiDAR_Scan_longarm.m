function [range_p, range_n] = LiDAR_Scan_longarm(map,pos,phi,...
    arm_len, show_p, show_n)
% Function: simulate a long-arm lidar in a 2D map (top view).
% Input:
%     map - array of 2D map (WxH), right & down are positive. 
%           The map defines distance unit.
%     pos - X, Y coordinates of the LiDAR in the map.
%     phi - array of long-arm orientaion (Nx1), deg.
%     arm_len - arm length.
%     show_p, show_n - time interval when showing laser (ms). 
%                  0 means no showing.
% Output:
%     range_p, range_n - array of range corresponding to phi (Nx1).
% Demo:
% 
% Writen by LIN, Jingyu (linjy02@hotmail.com), 20210121
%
if show_p>0 || show_n>0
    figure(100); imagesc(map); colormap(gray(256))
    hold on
end
N = length(phi);
d_e = 0.2;

% start scanning
x_pos = pos(1); y_pos = pos(2);
range_p = zeros(N,1);
range_n = zeros(N,1);
for i = 1:N
    x_phi = x_pos + arm_len*cosd(phi(i));
    y_phi = y_pos - arm_len*sind(phi(i)); % down is positive
    % front scannning
    [x_p,y_p,range_p(i)] = onMap_laser(x_phi,y_phi,phi(i) + 90,map,d_e);
    if show_p>0
        figure(100); plot([x_phi x_p],[y_phi y_p],'r-')
        pause(show_p/1000)
        drawnow
    end
    % back scanning
    [x_n,y_n,range_n(i)] = onMap_laser(x_phi,y_phi,phi(i) - 90,map,d_e);
    if show_n>0
        figure(100); plot([x_phi x_n],[y_phi y_n],'r-')
        pause(show_n/1000)
        drawnow
    end
end
if show_p>0 || show_n>0
    figure(100); hold off
end
