function [pnt_dens] = KDE_2D(pos, pntcloud, hG) 
% Function: estimate point density at a point by kernel density estimation.
% Input:
%     pos - x, y coordinates of the target point (m)
%     pntcloud - x, y coordinates of source points (K * 2)
%     hG - width of Gaussian (m)
% Output:
%     pnt_dens - density at the target point (points/m^2). 
% 
% Writen by LIN, Jingyu (linjy02@hotmail.com), 20210115
%
h_2 = hG*hG;
dx = pntcloud(:,1) - pos(1);
dy = pntcloud(:,2) - pos(2);
d2 = (dx.*dx+dy.*dy);
dens = exp(-d2/h_2);
pnt_dens = sum(dens(:))/hG/sqrt(pi);
