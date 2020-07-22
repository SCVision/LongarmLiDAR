function pnt_new = RotationZ(pnt_origin, alpha) 
% Function: general coordinates rotation  around Z-axis.
% Input:
%     pnt_origin - x, y, z coordinates of points (K * 3)
%     alpha - angle to rotate (deg)
% Output:
%     pnt_new - x, y, z coordinates of new points (K * 3)
% Demo: 
% alpha = 45;
% N = 10; % number of points in a line
% Xline = [(0:N)', zeros(N+1,1), zeros(N+1,1)];
% Yline = [zeros(N*5+1,1), (0:0.2:N)', zeros(N*5+1,1)];
% Zline = [zeros(N*10+1,1), zeros(N*10+1,1), (0:0.1:N)'];
% pnt_origin = [Xline;Yline;Zline];
% pnt_new = RotationZ(pnt_origin, alpha);
% figure(1); scatter3(pnt_origin(:,1),pnt_origin(:,2),pnt_origin(:,3),100,'.'); 
% xlabel('x'); ylabel('y'); zlabel('z');
% figure(2); scatter3(pnt_new(:,1),pnt_new(:,2),pnt_new(:,3),100,'.'); 
% xlabel('x'); ylabel('y'); zlabel('z');
% 
% Writen by LIN, Jingyu (linjy02@hotmail.com), 20200715
%
Rot = [cosd(alpha) -sind(alpha) 0; sind(alpha) cosd(alpha) 0; 0 0 1];
pnt_new = pnt_origin*Rot';
