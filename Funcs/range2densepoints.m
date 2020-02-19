function [pntCloud,map3D,X,Y,Z] = range2densepoints(range, angleV, angleH, R, Dtheta, ...
    box, boxStep, minRange, maxRange, maxDrange) 
% Function: generate dense point cloud in a box region from range mesh.
% Input:
%     range - range mesh data (H*V). 
%     angleV - vertical angles theta (V*1).
%     angleH - horizontal angles phi (H*1). 
%     R, Dtheta - calibration parameters 
%     box - definition of the box: [x1 y1 z1; x2 y2 z2] 
%     boxStep - point interval in the box
%     minRange, maxRange, maxDrange - valid range scope and delta range
% Output:
%     pnts - x, y, z coordinates of points ([X1 Y1 Z1; X2 Y2 Z2; ...])
%     map3D - x,y,z grid, 1=block, 0.5=unknown, 0=clear.
%     X,Y,Z - coordinates of map3D;
% Demo:
% 
% Writen by LIN, Jingyu (linjy02@hotmail.com), 20200215

% initializing
angleV = angleV + Dtheta;
x1 = round(box(1,1)/boxStep)*boxStep;
x2 = round(box(2,1)/boxStep)*boxStep; 
if x1<x2
    X = x1:boxStep:x2;
else
    X = x2:boxStep:x1;
end
y1 = round(box(1,2)/boxStep)*boxStep; 
y2 = round(box(2,2)/boxStep)*boxStep; 
if y1<y2
    Y = y1:boxStep:y2;
else
    Y = y2:boxStep:y1;
end
z1 = round(box(1,3)/boxStep)*boxStep; 
z2 = round(box(2,3)/boxStep)*boxStep;
if z1<z2
    Z = z1:boxStep:z2;
else
    Z = z2:boxStep:z1;
end

% main loop
thick = boxStep; % raduis of endpoint of laser beam
inner2 = minRange*minRange; 
R2 = R*R;
map3D = 0.5*ones(length(X),length(Y),length(Z));
pnts = zeros(numel(map3D), 3); % point cloud 
p_pnts = 0; % pointer to pnts
for i = 1:length(X)
    for j = 1:length(Y)
        D0 = X(i)*X(i) + Y(j)*Y(j) - R2;
        if D0<inner2 % inner region
            continue;
        end
        beta = atand(sqrt(D0)/R);
        alpha = atan2d(X(i),Y(j));
        phi1 = alpha - beta;
        phi2 = alpha + beta;
        for k = 1:length(Z)
            d = sqrt(D0 + Z(k)*Z(k)); % range of (x(i),y(j),z(k))          
            theta1 = acosd(Z(k)/d);
            theta2 = -theta1;
            
            % check if (x(i),y(j),z(k)) on a surface
            [r1, err1] = interpolate_mesh(theta1, phi1, range, angleV, angleH, minRange, maxRange, maxDrange);           
            [r2, err2] = interpolate_mesh(theta2, phi2, range, angleV, angleH, minRange, maxRange, maxDrange);           
            % four events: 1.error, 2.<d-boxStep, 3.inside, 4.>d+boxStep
            % [possible table]     r1: (err) (<) (in)  (>)
            % [? - unknown   ]    (err)  ?    ?    S    C
            % [S - surface   ] r2: (<)   ?    ?    S    C
            % [C - clear     ]     (in)  S    S    S    S
            %                      (>)   C    C    S    C
            if (r1 > d - thick && r1 < d + thick) || ...
                    (r2 > d - thick && r2 < d + thick)
                % surface point
                map3D(i,j,k) = 1;
                p_pnts = p_pnts + 1;
                pnts(p_pnts,:) = [X(i),Y(j),Z(k)];
            elseif r1 > d + thick && r2 > d - thick
                % clear point
                map3D(i,j,k) = 0;
            else
                % unknown region
            end
        end
    end
end
pntCloud = pnts(1:p_pnts,:);
