function pntcloud = range2points(range, angleV, angleH, R, Dtheta) 
% Function: convert range array to point sets.
% Input:
%     range - range data (H*V). 
%     angleV - vertical angles theta (V*1).
%     angleH - horizontal angles phi (H*1). 
%     R, Dtheta - calibration parameters 
% Output:
%     pntcloud - x, y, z coordinates of points ((H*V) * 3)
% Demo:
% R = 0.2; Dtheta = 0;
% [range, angleV, angleH, timestamp] = read_scandata('Scanned1.txt'); 
% ps = Rang2Points(range, angleV, angleH, R, Dtheta);
% figure(1); 
% scatter3(ps(:,1),ps(:,2),ps(:,3),1,'.'); xlabel('x'); ylabel('y'); zlabel('z'); 
% 
% Writen by LIN, Jingyu (linjy02@hotmail.com), 20200202
%
angleV = angleV(:)' + Dtheta; % row vector
SV = sind(angleV); % row vector
CV = cosd(angleV); % row vector
SH = sind(angleH); % row vector or column vector
CH = cosd(angleH); % row vector or column vector
[n,m] = size(range);
x = zeros(n,m);
y = zeros(n,m);
z = zeros(n,m);
for i = 1:n
    x(i,:) =  range(i,:).*SV*CH(i)+ R*SH(i);
    y(i,:) = -range(i,:).*SV*SH(i)+ R*CH(i);
    z(i,:) =  range(i,:).*CV;
end
pntcloud = [x(:), y(:), z(:)];