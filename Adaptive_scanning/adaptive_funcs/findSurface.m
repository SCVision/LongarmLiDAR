function [x_obj, z_obj, objDens, beta] = findSurface(pln_dens, xL, zL, Dphi_s0, Dtheta_s0) 
% Function: find object surface based on density.
% Input:
%     pnt_dens - density array (points/m^2) (M*N). 
%     xL - x coordinates on the plane (m) (N*1)
%     zL - z coordinates on the plane (m) (M*1)
%     Dphi_s0 - horizontal angular step of LIDAR (deg)
%     Dtheta_s0 - vertical angular step of LIDAR (deg)
% Output:
%     (x_obj, z_obj) - object coordinates on the scanning plane
%     objDens - density of points on object surface (points/m^2)
%     beta - slope angle of object surface (deg) 
% Test Demo:
% 
% Writen by LIN, Jingyu (linjy02@hotmail.com), 20200622
%
rat_null = 0.26; % null space when below this ratio of density to theoretical

xDens = sum(pln_dens,1);
[~, idx] = max(xDens);
x_obj = xL(idx); % distance from 2D LIDAR to object
[objDens, idx] = max(pln_dens(:,idx)); % density of object surface
z_obj = zL(idx);
% rm_2 = dm*dm + zm*zm;
Kdr = pi/180; % convert degree to radian
patch = x_obj*x_obj*Dphi_s0*Kdr*Dtheta_s0*Kdr; % 1/patch is theoretical density 
rat = objDens*patch; % objDens/(1/patch)
beta = 90;
if rat < rat_null % ratio is very small, maybe null space 
    x_obj = 0; z_obj = 0;
    return;
elseif rat >1 % that is impossible
    return;
else 
    beta = asind(rat); % slope angle
end
% e_m = 1/objDens/dm/(Dtheta_s0*Kdr);