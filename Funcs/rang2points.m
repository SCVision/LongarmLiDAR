function [X, Y, Z] = rang2points(range, angleV, angleH, R, Dtheta) 
% Function: convert range array to point sets.
% Input:
%     range - range data (H*V). 
%     angleV - vertical angles theta (V*1).
%     angleH - horizontal angles phi (H*1). 
%     R, Dtheta - calibration parameters 
% Output:
%     X, Y, Z - coordinates of points 
% Demo:
% R = 0.2; Dtheta = 0;
% [range, angleV, angleH, timestamp] = read_scandata('Scanned1.txt'); 
% [X, Y, Z] = Rang2Points(range, angleV, angleH, R, Dtheta);
% 
% Writen by LIN, Jingyu (linjy02@hotmail.com), 20200127
%

[n,m] = size(range);
%%%%%%%%%%%%%%%calculate coordinate%%%%%%%%%%%%%%%%%
x = zeros(n,m);
y = zeros(n,m);
z = zeros(n,m);

for i = 1:1:n
   for j = 1:1:m
       x(i,j) =  range(i,j)*sind(angleV(j)+Dtheta)*cosd(angleH(i))+ R*sind(angleH(i));
       y(i,j) = -range(i,j)*sind(angleV(j)+Dtheta)*sind(angleH(i))+ R*cosd(angleH(i));
       z(i,j) =  range(i,j)*cosd(angleV(j)+Dtheta); 
   end
end
X = reshape(x', n*m, 1);
Y = reshape(y', n*m, 1);
Z = reshape(z', n*m, 1);