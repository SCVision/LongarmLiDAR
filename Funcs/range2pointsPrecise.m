function pntcloud = range2pointsPrecise(range, angleV, angleH, ...
    R, Dphi, Dpsi, Dtheta) 
% Function: convert range array to point sets.
% Method: use z-x-y nautical angles 
% Input:
%     range - range data (H*V). 
%     angleV - vertical angles theta (V*1).
%     angleH - horizontal angles phi (H*1). 
%     R,Dtheta1, Dphi1, Dtheta2 - calibrated parameters 
% Output:
%     pntcloud - x, y, z coordinates of points ((H*V) * 3)
% Demo:
% R = 0.2; Dphi = 0; Dpsi = 0; Dtheta = 0;
% [range, angleV, angleH, timestamp] = read_scandata('Scanned1.txt'); 
% ps = range2pointsPrecise(range, angleV, angleH, R, Dphi, Dpsi, Dtheta);
% figure(1); 
% scatter3(ps(:,1),ps(:,2),ps(:,3),1,'.'); xlabel('x'); ylabel('y'); zlabel('z'); 
% 
% Writen by LIN, Jingyu (linjy02@hotmail.com), 20200716
%

% preprocessing
angleV = angleV(:)' + Dtheta; % first rotation
SV = sind(angleV); % row vector
CV = cosd(angleV); % row vector
% Yrot = [cosd(Dphi) 0 sind(Dphi); 0 1 0; -sind(Dphi) 0 cosd(Dphi)];
% Zrot = [cosd(Dpsi) -sind(Dpsi) 0; sind(Dpsi) cosd(Dpsi) 0; 0 0 1];
% Rot = Yrot*Zrot; % rotation from deviation angle
Zrot = [cosd(Dphi) -sind(Dphi) 0; sind(Dphi) cosd(Dphi) 0; 0 0 1];
Xrot = [1 0 0; 0 cosd(Dpsi) -sind(Dpsi); 0 sind(Dpsi) cosd(Dpsi)];
Rot = Zrot*Xrot; % rotation from deviation angle

% prepare for transform
[H,V] = size(range);
x = zeros(H,V);
y = zeros(H,V);
z = zeros(H,V);
zeroLines = zeros(1,V);
SH = sind(angleH); 
CH = cosd(angleH); 
for i = 1:H % for each scanning plane
    % polar coordinates to Cartesian coordinates
    xL_hat = range(i,:).*SV;
    zL_hat = range(i,:).*CV;
    
    % scanning plane to LIDAR coordinates
    x_tilde = Rot*[xL_hat;zeroLines;zL_hat]; % rotation from deviation angle
    x_tilde(2,:) = x_tilde(2,:) + R; % translation
    x(i,:) =  x_tilde(1,:)*CH(i)+ x_tilde(2,:)*SH(i);
    y(i,:) = -x_tilde(1,:)*SH(i)+ x_tilde(2,:)*CH(i);
    z(i,:) =  x_tilde(3,:);
end
pntcloud = [x(:), y(:), z(:)];
