function show_pointcloud(range, angleV, angleH, R, DTheta) 
% Function: show pointcloud with raw data scanned by our 3D Lidar.
%     range - range data (H*V). 
%     angleV - vertical angles (V lines).
%     angleH - horizontal angles (H lines). 
% Demo:
% [range, angleV, angleH, timestamp] = read_scandata('Scanned1.txt'); 
% show_pointcloud(range, angleV, angleH);
% 
% Writen by LI, Shuqing (jay_issaac@163.com), 20191115
%

[n,m] = size(range);
%%%%%%%%%%%%%%%calculate coordinate%%%%%%%%%%%%%%%%%
angleV = angleV+DTheta;
x = zeros(n,m);
y = zeros(n,m);
z = zeros(n,m);
for i = 1:1:n
   for j = 1:1:m
       x(i,j) = range(i,j)*sind(angleV(j))*cosd(angleH(i))+ R*sind(angleH(i));%
       y(i,j) = -range(i,j)*sind(angleV(j))*sind(angleH(i))+ R*cosd(angleH(i));%
       z(i,j) = range(i,j)*cosd(angleV(j));       
   end
end
%%%%%%%%%%%%%%%draw pointcloud%%%%%%%%%%%%%%%%%
X = reshape(x.',n*m,1);
Y = reshape(y.',n*m,1);
Z = reshape(z.',n*m,1);
% figure();
scatter3(X,Y,Z,1,'.');  
xlabel('x');
ylabel('y');
zlabel('z');
% title('pointcloud display');