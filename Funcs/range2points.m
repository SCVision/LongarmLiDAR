function points = range2points(range, angleV, angleH, R, Dtheta) 
% Function: convert range array to point sets.
% Input:
%     range - range data (H*V). 
%     angleV - vertical angles theta (V*1).
%     angleH - horizontal angles phi (H*1). 
%     R, Dtheta - calibration parameters 
% Output:
%     points - X, Y, Z coordinates of points ((H*V) * 3)
% Demo:
% R = 0.2; Dtheta = 0;
% [range, angleV, angleH, timestamp] = read_scandata('Scanned1.txt'); 
% ps = Rang2Points(range, angleV, angleH, R, Dtheta);
% figure(1); 
% scatter3(ps(:,1),ps(:,2),ps(:,3),1,'.'); xlabel('x'); ylabel('y'); zlabel('z'); 
% 
% Writen by LIN, Jingyu (linjy02@hotmail.com), 20200202
%

SV = sind(angleV(:)'+Dtheta);
CV = cosd(angleV(:)'+Dtheta);
SH = sind(angleH);
CH = cosd(angleH);
[n,m] = size(range);
x = zeros(n,m);
y = zeros(n,m);
z = zeros(n,m);
for i = 1:1:n
    x(i,:) =  range(i,:).*SV*CH(i)+ R*SH(i);
    y(i,:) = -range(i,:).*SV*SH(i)+ R*CH(i);
    z(i,:) =  range(i,:).*CV;
end
% X = reshape(x', n*m, 1);
% Y = reshape(y', n*m, 1);
% Z = reshape(z', n*m, 1);
points = [x(:), y(:), z(:)];