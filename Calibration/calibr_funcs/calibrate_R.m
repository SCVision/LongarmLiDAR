function [R,error]=calibrate_R(d1, theta1, phi1, d2, theta2, phi2, Dtheta, iter) 
% Function: Calibrate R. 
% Setup: The LIDAR is horizontal and a calibration board is vertical.
%    Scan by left side and right side of the 2D LIDAR on the same plane
%    of Ax+By+Cz=1, respectively. 
% Input:
%    d1, d2 - range data (H*V). 
%    theta1, theta2 - vertical angles (V*1). They have different signs.
%    phi1, phi2 - horizontal angles (H*1).
%    Dtheta - calibrated delta-theta
%    iter - iteration times. 
% Output:
%    R - calibrated parameters. 
%    error - cabration error. 
% 
% Writen by LIN Jingyu (linjy02@hotmail.com), 20200129
%

% prepare constants
theta1 = theta1 + Dtheta; theta2 = theta2 + Dtheta;

[n,m] = size(d1);
K11=zeros(n,m);
K12=zeros(n,m);
K13=zeros(n,m);
K14=zeros(n,m);
K15=zeros(n,m);
for i = 1:n
    for j = 1:m
        K11(i,j) = d1(i,j)*sind(theta1(j))*cosd(phi1(i));
        K12(i,j) = d1(i,j)*sind(theta1(j))*sind(phi1(i));
        K13(i,j) = sind(phi1(i));
        K14(i,j) = cosd(phi1(i));
        K15(i,j) = d1(i,j)*cosd(theta1(j));
    end
end

[n,m] = size(d2);
K21=zeros(n,m);
K22=zeros(n,m);
K23=zeros(n,m);
K24=zeros(n,m);
K25=zeros(n,m);
for i = 1:n
    for j = 1:m
        K21(i,j) = d2(i,j)*sind(theta2(j))*cosd(phi2(i));
        K22(i,j) = d2(i,j)*sind(theta2(j))*sind(phi2(i));
        K23(i,j) = sind(phi2(i));
        K24(i,j) = cosd(phi2(i));
        K25(i,j) = d2(i,j)*cosd(theta2(j));
    end
end

% iteration
R = 0;
min_error = 1e6;
error = zeros(iter,1);
err_p = 1;
lr = 0.001;
for ii=1:iter
    R = R + lr;

    % calculate A, B, C
    Kx = [K11(:)+R*K13(:), -K12(:)+R*K14(:), K15(:)];
    X1 = Kx \ ones(size(d1(:))); % X = [A;B;C]
    Kx = [K21(:)+R*K23(:), -K22(:)+R*K24(:), K25(:)];
    X2 = Kx \ ones(size(d2(:))); % X = [A;B;C]
    
    % record error
    a = log(sum((X1 - X2).^2));
    if a > min_error 
        % search back if error increase
        lr = -lr*0.618;
        if abs(lr) <1e-10
            break
        end
    end
    min_error = a;
    error(err_p) = a;
    err_p = err_p + 1;
end
error = error(1:err_p-1);
% fprintf('R = %f, error = %f\n', R, error(end));