% Show data in 'anti_occlusion_demo'
path(path,'..\Funcs')

%% read data
dir = 'anti_occlusion_demo';
dataset = 4;
switch dataset
    case 1
        subdir = '1r-4-350-1-1050-2-350-3-1050';
        fn = 'batchScanned20200705091246.txt';
    case 2
        subdir = '1r-4-350-1-1050-2-350-3-1050';
        fn = 'batchScanned20200705092021.txt';
    case 3
        subdir = '1r-4-350-1-1050-2-350-3-1050';
        fn = 'batchScanned20200705092809.txt';
    otherwise
        subdir = '1r-4-350-1-1050-2-350-3-1050-basket';
        fn = 'batchScanned20200705094609.txt';
end
[range, angleV, angleH, timestamp] = read_scandata([dir,'\',subdir,'\',fn]); 
range = replace_outlier(range,0.1, 10);

%% show point cloud
R=0.195000; Dphi=0.0000; Dpsi=0.4000; Dtheta=-2.2000;
ps = range2pointsPrecise_ts(range, angleV, angleH, R,Dphi,Dpsi,Dtheta);
figure(1); 
scatter3(ps(:,1),ps(:,2),ps(:,3),1,'.');
az = 190; el = 60; view(az,el)
xlim([-2,2]); ylim([-1.5 1.5]); 
zlim([-0.16,0.5])
xlabel('x'); ylabel('y'); zlabel('z'); 

