function [pln_dens] = KDEonVplane(xyp, zp, pntcloud, hG) 
% Function: kernel density estimation on a vertival plane.
% Input: 
%     xyp - x, y coordinates on the plane (m) (N*2)
%     zp - z coordinates on the plane (m) (M*1)
%     pntcloud - x, y, z coordinates of source points (K * 3)
%     h - width of Gaussian (m)
% Output:
%     pln_dens - density (points/m^2) (M*N). 
% Test Demo: none
% 
% Writen by LIN, Jingyu (linjy02@hotmail.com), 20200622
%
M = length(zp); N = length(xyp);
pln_dens = zeros(M, N);
for i=1:M
    for j=1:N
        pos = [xyp(j,1),xyp(j,2),zp(i)];
        pln_dens(i,j) = KDE(pos, pntcloud, hG);
    end
end
