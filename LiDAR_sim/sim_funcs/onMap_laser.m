function [X,Y,dist] = onMap_laser(Xc,Yc,theta,map,d_e,range_out)
% Function: simulate a beam of laser in a 2D map.
% Input:
%     Xc,Yc - starting point of laser in the map.
%     theta - direction of laser, deg.
%     map - array of 2D map (WxH), right & down are positive. 
%           The map defines distance unit.
%     d_e - precision of range.
%     range_out - range for out of range
% Output:
%     X,Y - position of the object being found.
%     dist - distance of the object being found, large value for out of map
% 
% Writen by LIN, Jingyu (linjy02@hotmail.com), 20210121
%

% map information
map = double(map);
[H, W] = size(map);
max_rang = sqrt(H*H+W*W);

if nargin < 5, d_e = 0.2; end  
if nargin < 6, range_out = 1e8; end  

% emit beam
for t = 1:d_e:max_rang
    X = Xc + t*cosd(theta);
    Y = Yc - t*sind(theta); % down is positive
    % check border
    if Y<1 || Y>H || X<1 || X>W
        % reach border
        dist = range_out;
        X = Xc + range_out*cosd(theta);
        Y = Yc - range_out*sind(theta); % down is positive
        break;
    end
    if map(round(Y),round(X)) <= 0
        % reach objects
        dist = t;
        break;
    end
end
