function [val, err] = interpolate_mesh(x, y, mesh_data, Xs, Ys, ...
    minVal, maxVal, dVal)
% Function: get a value bilinear-interpolated from a mesh.
% Input:
%     (x,y) - coordinates in mesh_data, not column and row position. 
%     mesh_data - data (Y*X). 
%     Xs - x coordinats, sorted (X*1).
%     Ys - y coordinats, sorted (Y*1). 
%     maxVal, minVal, dVal - valid value and delta value in mesh_data
% Output:
%     val - interpolated value. -1 if out of scope.
%     err - 'true' means out of scope.
%
% 
% Writen by LIN, Jingyu (linjy02@hotmail.com), 20200215
%

val = -1; err = true;
% Xs and Ys must be sorted from small to great. 
if Xs(1) > Xs(end)
    Xs = -Xs;
    x = -x;
end
if Ys(1) > Ys(end)
    Ys = -Ys;
    y = -y;
end
% search in mesh
p_x = find_position(Xs, x);
p_y = find_position(Ys, y);
if p_x < 0 || p_y < 0
    return; % err = true;
end

% calculate weight
errorfound = 0;
dx1 = x - Xs(p_x);
dx2 = Xs(p_x+1) - x;
dy1 = y - Ys(p_y);
dy2 = Ys(p_y+1) - y;
r11 = mesh_data(p_y, p_x); 
r12 = mesh_data(p_y, p_x+1); 
if abs(r11 - r12) > dVal
    return; % err = true;
end
if r11 > minVal && r11 < maxVal && r12 > minVal && r12 < maxVal
    r01 = (r11*dx2+r12*dx1)/(dx1+dx2);
else
    errorfound = 1;
    if r11 > minVal && r11 < maxVal
        r01 = r11;
    elseif r12 > minVal && r12 < maxVal
        r01 = r12;
    else
        return; % err = true;
    end
end
r21 = mesh_data(p_y+1, p_x); 
r22 = mesh_data(p_y+1, p_x+1); 
if abs(r21 - r22) > dVal
    return; % err = true;
end
if r21 > minVal && r21 < maxVal && r22 > minVal && r22 < maxVal
    r02 = (r21*dx2+r22*dx1)/(dx1+dx2);
else
    if errorfound > 0
        return; % err = true;
    elseif r21 > minVal && r21 < maxVal
        r02 = r21;
    elseif r22 > minVal && r22 < maxVal
        r02 = r22;
    else
        return; % err = true;
    end
end
if abs(r01 - r02) > dVal
    return; % err = true;
end
val = (r01*dy2+r02*dy1)/(dy1+dy2);
err = false;
return

function xPos = find_position(Xcoordinates, xVal)
% Function: find the postion of 'xVal' in sorted 'coordinates'
% Output: xPos satisfies xVal in [Xcoordinates(xPos), Xcoordinates(xPos+1)).
% Note: Xcoordinates must be sorted from small to great.
if xVal < Xcoordinates(1) || xVal > Xcoordinates(end)
    xPos = -1; % not found
    return;
end
N = length(Xcoordinates);
p1 = 1;
p2 = N;
dataLen = p2 - p1 + 1;
xPos = floor(dataLen / 2);
while dataLen >= 3 && xPos < N
    a = Xcoordinates(xPos);
    if xVal >= a
        if xVal < Xcoordinates(xPos+1)
            break; % position found
        end
        p1 = xPos;
    else % theta < a
        p2 = xPos;
    end
    dataLen = p2 - p1 + 1;
    xPos = p2 - floor(dataLen / 2);
end
return