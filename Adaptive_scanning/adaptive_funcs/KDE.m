function [pnt_dens] = KDE(pos, pntcloud, hG) 
% Function: kernel density estimation.
% Input:
%     pos - x, y, z coordinates of positions (m)
%     pntcloud - x, y, z coordinates of source points (K * 3)
%     hG - width of Gaussian (m)
% Output:
%     pnt_dens - density at each position (points/m^2). 
% Test Demo:
% x = 0:0.1:10; y = 0:20; % density = 10/squre unit
% h = 1;
% [X,Y] = meshgrid(x,y);
% pntcloud = [X(:) Y(:) 0*ones(length(X(:)),1)];
% figure(1); 
% scatter3(pntcloud(:,1),pntcloud(:,2),pntcloud(:,3),20,'.'); 
% xlabel('x'); ylabel('y'); zlabel('z'); 
% % calculate density
% pnt_dens = zeros(length(y), length(x));
% z = 0;
% for i=1:length(y)
%     for j=1:length(x)
%         pos = [x(j),y(i),z];
%         pnt_dens(i,j) = KDE(pos, pntcloud, h);
%     end
% end
% figure(2); mesh(X,Y,pnt_dens); 
% xlabel('x'); ylabel('y')
% 
% Writen by LIN, Jingyu (linjy02@hotmail.com), 20200622
%
if length(pos) ~= 3
    pnt_dens = 0;
    return
end
h_2 = hG(1)*hG(1);
dx = pntcloud(:,1) - pos(1);
dy = pntcloud(:,2) - pos(2);
dz = pntcloud(:,3) - pos(3);
d2 = (dx.*dx+dy.*dy+dz.*dz);
dens = exp(-d2/h_2);
pnt_dens = sum(dens(:))/h_2/pi;