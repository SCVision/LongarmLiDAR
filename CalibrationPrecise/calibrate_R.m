function [R, err] = calibrate_R(d1,theta1,phi1, d2,theta2,phi2,...
    R, Dphi, Dpsi, Dtheta, iter, Rstep, msgOn) 
% Function: calibrate long-arm LIDAR parameters.
% Method: use z-x-y nautical angles 
% Method: calibrate R on a wall
% Input:
%     d1, d2 - front and back range data (H*V). 
%     theta1, theta2 - vertical angles theta (V*1).
%     phi1, phi2 - horizontal angles phi (H*1). 
%     R, Dphi, Dpsi, Dtheta - initial calibrated parameters
%     iter - iteration times
%     Rstep - step for searching R
%     msgOn - true to display progress
% Output:
%     R, err - calibrated parameters and error
%     err - calibration error in plane fitting
% Demo: none
% 
% Writen by LIN, Jingyu (linjy02@hotmail.com), 20200718
% 
[H1,V1] = size(d1); % front data
[H2,V2] = size(d2); % back data
K = H1*V1 + H2*V2; % number of points
err = zeros(iter,1); % error for 2 calibration boards

% calibrate R
min_error = 1e6;
minRstep = Rstep*1e-4; % meter
for ii = 1:iter
    R = R + Rstep;
   
    % calculate c = [A;B;C]
    ps1 = range2pointsPrecise(d1,theta1,phi1,R,Dphi,Dpsi, Dtheta);
    ps2 = range2pointsPrecise(d2,theta2,phi2,R,Dphi,Dpsi, Dtheta);
    X = [ps1;ps2];
    c = X\ones(K,1);
    err_plane = X*c - 1;
    err(ii) = err_plane'*err_plane/2/K;
    if err(ii) > min_error 
        % search back if error increase
        Rstep = -Rstep*0.618;
        if abs(Rstep) < minRstep
            break
        end
    end
    min_error = err(ii);
    
    if msgOn
        figure(101); plot(ii, log(err(ii)), 'r*');
        title(' R - log error')
        xlim([0,iter]); hold on
    end
end
figure(101); hold off
if msgOn
    fprintf('R=%.6f; Dphi=%.4f; Dpsi=%.4f; Dtheta=%.4f;\n',R,Dphi,Dpsi,Dtheta)
end
