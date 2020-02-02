function [Dtheta,b]=calibrate_Dtheta(d, theta) 
% Function: Calibrate theta when LIDAR system is parallel to calibration board,
%           therefore scanning points in the line zL=k*xL+b, k=tan(Dtheta)
% Input:
%     d - range data (H*V). 
%     theta - vertical angles (V lines).
% Output:
%     Dtheta - calibrated angle. 
%     b - parameter of line. 
% Demo:
% [range, angleV, angleH, timestamp] = read_scandata('Scanned1.txt'); 
% d=range(1:1,:);
% theta=angleV(:);
% [Dtheta,b]=calibrate_theta(d, theta);
% 
% Writen by LI, Shuqing (jay_issaac@163.com), 20191115
%

A=[d.*sind(theta)';ones(1,size(theta,1))]';
B=(d.*cosd(theta)')';
result=A\B;
k=result(1);
Dtheta=atand(k);
b=result(2);
fprintf('Dtheta = %f, b = %f\n',Dtheta, b);
