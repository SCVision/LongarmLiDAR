function [Dtheta, err] = calibrate_Dtheta(d_ceil, theta_ceil, phi_ceil, ...
    R, Dphi, Dpsi, Dtheta, iter, Dtheta_step, msgOn) 
% Function: calibrate long-arm LIDAR parameters.
% Method: use z-x-y nautical angles 
% Method: calibrate Dtheta on a ceiling
% Input:
%     d_ceil - ceiling range data (H*V). 
%     theta_ceil - vertical angles theta (V*1).
%     phi_ceil - horizontal angles phi (H*1). 
%     R, Dphi, Dpsi, Dtheta - initial calibrated parameters
%     iter - iteration times
%     Dtheta_step - step for searching Dtheta
%     msgOn - true to display progress
% Output:
%     Dtheta, err - calibrated parameters and error
%     err - calibration error in plane fitting
% Demo: none
% 
% Writen by LIN, Jingyu (linjy02@hotmail.com), 20200719
% 
[H,V] = size(d_ceil);
K = H*V; % number of points
err = zeros(iter,1); % error for 2 calibration boards

% calibrate Dtheta
min_error = 1e6;
minDthetastep = Dtheta_step*1e-3; % deg
for ii = 1:iter
    Dtheta = Dtheta + Dtheta_step;
    
    % calculate c = [A;B;C]
    X = range2pointsPrecise(d_ceil,theta_ceil,phi_ceil,R,Dphi,Dpsi, Dtheta);
    c = X\ones(K,1);
    err_plane = X*c - 1;
    err(ii) = err_plane'*err_plane/2/K;
    if err(ii) > min_error 
        % search back if error increase
        Dtheta_step = -Dtheta_step*0.618;
        if abs(Dtheta_step) < minDthetastep
            break
        end
    end
    min_error = err(ii);
    
    if msgOn
        figure(100); plot(ii, log(err(ii)), 'r*');
        title(' Dtheta - log error')
        xlim([0,iter]); hold on
    end
end
figure(100); hold off
if msgOn
    fprintf('First calibration: Dtheta=%.4f;\n',Dtheta)
end
