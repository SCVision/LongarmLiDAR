function pntcloud = range2pointsPrecise_ts(rData, angleV, angleH, ...
    R, Dphi, Dpsi, Dtheta) 
% Function: convert range array to point sets.
% Method: use z-x-y nautical angles 
% Method: involve temporal factor 
% Input:
%     rData - range data (H*V). 
%     angleV - vertical angles theta (V*1).
%     angleH - horizontal angles phi (H*1). 
%     R,Dtheta1, Dphi1, Dtheta2 - calibrated parameters 
% Output:
%     pntcloud - x, y, z coordinates of points ((H*V) * 3)
% Demo:
% R = 0.2; Dphi = 0; Dpsi = 0; Dtheta = 0;
% [rData, angleV, angleH, timestamp] = read_scandata('Scanned1.txt'); 
% ps = range2pointsPrecise(rData, angleV, angleH, R, Dphi, Dpsi, Dtheta);
% figure(1); 
% scatter3(ps(:,1),ps(:,2),ps(:,3),1,'.'); xlabel('x'); ylabel('y'); zlabel('z'); 
% 
% Writen by LIN, Jingyu (linjy02@hotmail.com), 20200722
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
[H,V] = size(rData);
x = zeros(H,V);
y = zeros(H,V);
z = zeros(H,V);
zeroLines = zeros(1,V);
Dphi_s = diff(angleH(:)) / V; % angular step of each point 
Dphi_s = [Dphi_s;Dphi_s(end)]; 
% SH = sind(angleH); 
% CH = cosd(angleH); 
for i = 1:H % for each scanning plane
    phi = angleH(i) + ((1:V)-1)*Dphi_s(i); % platform positions
    SH = sind(phi); CH = cosd(phi);
    
    % polar coordinates to Cartesian coordinates
    xL_hat = rData(i,:).*SV;
    zL_hat = rData(i,:).*CV;
    
    % scanning plane to LIDAR coordinates
    x_tilde = Rot*[xL_hat;zeroLines;zL_hat]; % rotation from deviation angle
    x_tilde(2,:) = x_tilde(2,:) + R; % translation
    x(i,:) =  x_tilde(1,:).*CH+ x_tilde(2,:).*SH;
    y(i,:) = -x_tilde(1,:).*SH+ x_tilde(2,:).*CH;
    z(i,:) =  x_tilde(3,:);
end
pntcloud = [x(:), y(:), z(:)];
