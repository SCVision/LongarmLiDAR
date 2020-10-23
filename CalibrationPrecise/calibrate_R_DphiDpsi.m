function [R, Dphi, Dpsi, err] = calibrate_R_DphiDpsi(d1,theta1,phi1,d2,theta2,phi2, ...
    R, Dphi, Dpsi, Dtheta, iter, lr, msgOn) 
% Function: calibrate long-arm LIDAR parameters.
% Method: use z-x-y nautical angles 
% Method: simultaneously calibrate R, Dphi, Dpsi on a wall by gradient descent.
% Input:
%     d1, d2 - front and back range data (H*V). 
%     theta1, theta2 - vertical angles theta (V*1).
%     phi1, phi2 - horizontal angles phi (H*1). 
%     R, Dphi, Dpsi, Dtheta - initial calibrated parameters
%     iter - iteration times
%     lr - learning rate (deg)
%     msgOn - true to display progress
% Output:
%     R, Dphi, Dpsi - calibrated parameters and error
%     err - calibration error in plane fitting
% Demo: none
% 
% Writen by LIN, Jingyu (linjy02@hotmail.com), 20200720
% 
N_var = 3;

% for front data
[H1,V1] = size(d1);
ps1 = zeros(H1*V1,3); % holding points
dx1_alpha = zeros(H1*V1,3*N_var); % holding gradients

% for back data
[H2,V2] = size(d2);
ps2 = zeros(H2*V2,3);
dx2_alpha = zeros(H2*V2,3*N_var); % holding gradients

% initialize
K = H1*V1 + H2*V2; % number of points
err = zeros(iter,1); % error for 2 calibration boards

% calibrate R, Dphi, Dpsi
angleV = theta1(:)' + Dtheta; % first rotation
SV1 = sind(angleV);
CV1 = cosd(angleV);
angleV = theta2(:)' + Dtheta; % first rotation
SV2 = sind(angleV);
CV2 = cosd(angleV);
lambda = lr; % learning rate (deg)
minLambda = lr*1e-4;
for ii = 1:iter
    Xrot = [1 0 0; 0 cosd(Dpsi) -sind(Dpsi); 0 sind(Dpsi) cosd(Dpsi)];
    Zrot = [cosd(Dphi) -sind(Dphi) 0; sind(Dphi) cosd(Dphi) 0; 0 0 1];
    Rot = Zrot*Xrot; % rotation from deviation angle

    % calculate points and gradients - front
    zeroLine = zeros(1,V1); % zero lines
    for i = 1:H1 % for each scanning plane
        % scanning plane: polar coordinates to Cartesian coordinates
        xL_hat = d1(i,:).*SV1;
        zL_hat = d1(i,:).*CV1;
        x_tilde = Rot*[xL_hat;zeroLine;zL_hat]; % xL_hat to x_tilde

        % gradients
        dx_R = [sind(phi1(i)), cosd(phi1(i)), 0];
        dx_tilde_dphi = [-x_tilde(2,:);x_tilde(1,:);zeroLine];
        dx_dphi = RotationZ(dx_tilde_dphi', -phi1(i));
        dx_tilde_dpsi = Rot*[zeroLine;-zL_hat;zeroLine];
        dx_dpsi = RotationZ(dx_tilde_dpsi', -phi1(i));
        dx1_alpha((1:V1)+V1*(i-1),1:3) = repmat(dx_R,[V1 1]);
        dx1_alpha((1:V1)+V1*(i-1),4:6) = dx_dphi;
        dx1_alpha((1:V1)+V1*(i-1),7:9) = dx_dpsi;
        
        % x_tilde to x_world
        x_tilde(2,:) = x_tilde(2,:) + R; % translation
        x_world = RotationZ(x_tilde', -phi1(i)); 
        ps1((1:V1)+V1*(i-1),:) = x_world;
    end

    % calculate points and gradients - back
    zeroLine = zeros(1,V2); % zero lines
    for i = 1:H2 % for each scanning plane
        % scanning plane: polar coordinates to Cartesian coordinates
        xL_hat = d2(i,:).*SV2;
        zL_hat = d2(i,:).*CV2;
        x_tilde = Rot*[xL_hat;zeroLine;zL_hat]; % xL_hat to x_tilde
        
        % gradients
        dx_R = [sind(phi2(i)), cosd(phi2(i)), 0];
        dx_tilde_dphi = [-x_tilde(2,:);x_tilde(1,:);zeroLine];
        dx_dphi = RotationZ(dx_tilde_dphi', -phi2(i));
        dx_tilde_dpsi = Rot*[zeroLine;-zL_hat;zeroLine];
        dx_dpsi = RotationZ(dx_tilde_dpsi', -phi2(i));
        dx2_alpha((1:V2)+V2*(i-1),1:3) = repmat(dx_R,[V2 1]);
        dx2_alpha((1:V2)+V2*(i-1),4:6) = dx_dphi;
        dx2_alpha((1:V2)+V2*(i-1),7:9) = dx_dpsi;
        
        % x_tilde to x_world
        x_tilde(2,:) = x_tilde(2,:) + R; % translation
        x_world = RotationZ(x_tilde', -phi2(i)); 
        ps2((1:V2)+V2*(i-1),:) = x_world;
    end
    
    % calculate c = [A;B;C]
    X = [ps1;ps2];
    c = X\ones(K,1);
    err_plane = X*c - 1;
    err(ii) = err_plane'*err_plane/2/K;
    if ii>1 && err(ii)>err(ii-1)
        lambda = lambda*0.618;
        if lambda < minLambda
            break;
        end
    end
    if msgOn
        figure(101); plot(ii, log(err(ii)), 'r*');
        title('R, Dphi, Dpsi - log error')
        xlim([0,iter]); hold on
    end
    
    % gradient of error
    dX = [dx1_alpha;dx2_alpha];
    derr_R = sum(dX(:,1:3)*c.*err_plane)/K;
    derr_Dphi = sum(dX(:,4:6)*c.*err_plane)/K;
    derr_Dpsi = sum(dX(:,7:9)*c.*err_plane)/K;
    mag = sqrt(derr_R*derr_R+derr_Dphi*derr_Dphi+derr_Dpsi*derr_Dpsi);
    
    % update parameters
    R = R - lambda/mag*derr_R;
    Dphi = Dphi - lambda/mag*derr_Dphi;
    Dpsi = Dpsi - lambda/mag*derr_Dpsi;
end
figure(101); hold off
if msgOn
    fprintf('R=%.6f; Dphi=%.4f; Dpsi=%.4f; Dtheta=%.4f;\n',R,Dphi, Dpsi, Dtheta)
end
