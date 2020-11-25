% Show data in 'Data'
path(path,'..\Funcs')

%% select data 
datano = 1; % select data 1~8
switch(datano)
    case 1
        subdir = '0726';
        fn = 'batchScanned20201119190726.txt';
        xl=[-3,3]; yl=[-3 3]; zl=[-0.3,3];
    case 2
        subdir = '0748';
        fn = 'batchScanned20201119210748.txt';
        xl=[-3,3]; yl=[-0 3]; zl=[-0.3,3];
    case 3
        subdir = '1727';
        fn = 'batchScanned20201119191727.txt';
        xl=[-3,3]; yl=[-3 3]; zl=[-0.3,3];
    case 4
        subdir = '3635';
        fn = 'batchScanned20201119203635.txt';
        xl=[-3,3]; yl=[-0 3]; zl=[-0.3,3];
    case 5
        subdir = '4429';
        fn = 'batchScanned20201119184429.txt';
        xl=[-3,3]; yl=[-0 3]; zl=[-0.3,3];
    case 6
        subdir = '4731';
        fn = 'batchScanned20201119204731.txt';
        xl=[-3,3]; yl=[-0 3]; zl=[-0.3,3];    
    case 7
        subdir = '5517';
        fn = 'batchScanned20201119185517.txt';
        xl=[-3,3]; yl=[-0 3]; zl=[-0.3,3]; 
    case 8
        subdir = '5736';
        fn = 'batchScanned20201119205736.txt';
        xl=[-3,3]; yl=[-0 3]; zl=[-0.3,3]; 

    otherwise % put new data here
        dataname = 'batchScanned20200703172752';
end

% read data
% x1 = -2.5; x2 = 1;
% y1 = -1; y2 = 1.8;
% z1 = -0.3; z2 = 0.9;

[range, angleV, angleH, timestamp] = read_scandata([subdir,'\',fn]); 
%range = replace_outlier(range,0.1, 10);

%% show point cloud
R=0.1919; Dphi=0.0000; Dpsi=0.4000; Dtheta=-1.74;
ps = range2pointsPrecise(range, angleV, angleH, R,Dphi,Dpsi,Dtheta);
figure(1); 
scatter3(ps(:,1),ps(:,2),ps(:,3),1,'.');
az = 10; el = 20; view(az,el)
% xlim([-2,0.3]); ylim([-0.4 2]); 
% zlim([-0.16,0.2])
% xlim(xl); ylim(yl); zlim(zl)
xlabel('x'); ylabel('y'); zlabel('z'); 

