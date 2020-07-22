% Show data in 'anti_occlusion_test'
path(path,'..\Funcs')
dir = 'anti_occlusion_test';

%% read data
subdir = '1r-9-350-3-1050-8-350-10-1050';
fn = 'batchScanned20200621221312.txt';
x1 = -2.5; x2 = 1;
y1 = -1; y2 = 1.8;
z1 = -0.3; z2 = 0.9;

[range, angleV, angleH, timestamp] = read_scandata([dir,'\',subdir,'\',fn]); 
range = replace_outlier(range,0.1, 10);

%% show point cloud
if subdir(1)=='0'
    R=0.0387; Dphi=0.0000; Dpsi=0.4000; Dtheta=-1.84;
else
    R=0.1919; Dphi=0.0000; Dpsi=0.4000; Dtheta=-1.74;
end
ps = range2pointsPrecise_ts(range, angleV, angleH, R,Dphi,Dpsi,Dtheta);
figure(1); 
scatter3(ps(:,1),ps(:,2),ps(:,3),1,'.');
az = 30; el = 60; view(az,el)
xlim([-2,0.3]); ylim([-0.4 2]); 
zlim([-0.16,0.2])
xlabel('x'); ylabel('y'); zlabel('z'); 

