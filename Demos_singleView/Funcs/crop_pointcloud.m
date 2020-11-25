function [pointsin, pointsout] = crop_pointcloud(points, box) 
% Function: segment points inside and outside a box.
% Input:
%     points - X, Y, Z coordinates of points (n * 3).
%     box - definition of the box: [x1 y1 z1; x2 y2 z2] 
% Output:
%     pointsin - X, Y, Z coordinates of points in the box (n_in * 3)
%     pointsout - X, Y, Z coordinates of points out of the box ((n-n_in)*3)
% Demo:
% R = 0.2; Dtheta = 0;
% [range, angleV, angleH, timestamp] = read_scandata('Scanned1.txt'); 
% ps = Rang2Points(range, angleV, angleH, R, Dtheta);
% x1 = -5; x2 = 0;
% y1 = -5; y2 = 2;
% z1 = -5; z2 = 1;
% [psin, psout] = crop_pointcloud(ps, [x1, y1, z1; x2, y2, z2]) 
% figure(1); 
% scatter3(psin(:,1),psin(:,2),psin(:,3),1,'.'); xlabel('x'); ylabel('y'); zlabel('z'); 
% 
% Writen by LIN, Jingyu (linjy02@hotmail.com), 20200202
%
x1 = box(1,1); x2 = box(2,1); 
y1 = box(1,2); y2 = box(2,2); 
z1 = box(1,3); z2 = box(2,3);
n = length(points);
pointBin = zeros(n,3);
n_in = 0; n_out = 0;
for i = 1:n
    % check if the point in the box
    if points(i,1)>=x1 && points(i,1)<=x2 && ...
        points(i,2)>=y1 && points(i,2)<=y2 && ...
        points(i,3)>=z1 && points(i,3)<=z2
        n_in = n_in + 1;
        pointBin(n_in,1) = points(i,1); 
        pointBin(n_in,2) = points(i,2);
        pointBin(n_in,3) = points(i,3);
    else
        pointBin(n-n_out,1) = points(i,1); 
        pointBin(n-n_out,2) = points(i,2);
        pointBin(n-n_out,3) = points(i,3);
        n_out = n_out + 1;
    end
end
pointsin = pointBin(1:n_in,:);
pointsout = pointBin(n:-1:n-n_out+1,:);