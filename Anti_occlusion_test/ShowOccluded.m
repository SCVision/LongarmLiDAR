% Show data in 'anti_occlusion_test'
path(path,'..\Funcs')
dir = 'anti_occlusion_test';

%% select data 
datano = 1; % select data 1~36
switch(datano)
    case 1
        subdir = '0r-9-350-1-700-8-350-3-700';
        fn = 'batchScanned20200620214518.txt';
        xl=[-0.2,1.5]; yl=[-1.5 0.4]; zl=[-0.16,1.5];
    case 2
        subdir = '0r-9-350-1-1050-8-350-3-1050';
        fn = 'batchScanned20200620215518.txt';
        xl=[-0.2,1.5]; yl=[-1.5 0.4]; zl=[-0.16,1.5];
    case 3
        subdir = '0r-9-350-3-700-8-350-10-700';
        fn = 'batchScanned20200621233005.txt';
        xl=[-0.2,1.5]; yl=[-1.5 0.4]; zl=[-0.16,1.5];
    case 4
        subdir = '0r-9-350-3-1050-8-350-10-1050';
        fn = 'batchScanned20200621232604.txt';
        xl=[-0.2,1.5]; yl=[-1.5 0.4]; zl=[-0.16,1.5];
    case 5
        subdir = '0r-9-520-1-1050-8-520-3-1050';
        fn = 'batchScanned20200620221130.txt';
        xl=[-0.2,1.5]; yl=[-1.5 0.4]; zl=[-0.16,1.5];
    case 6
        subdir = '0r-9-520-1-1400';
        fn = 'batchScanned20200621233525.txt';
        xl=[-0.2,1.5]; yl=[-1.5 0.4]; zl=[-0.16,1.5];    
    case 7
        subdir = '0r-9-520-3-1050-8-520-10-1050';
        fn = 'batchScanned20200621230732.txt';
        xl=[-0.2,1.5]; yl=[-1.5 0.4]; zl=[-0.16,1.5];  
    case 8
        subdir = '0r-9-520-3-1400-8-520-10-1400';
        fn = 'batchScanned20200621231221.txt';
        xl=[-0.2,1.5]; yl=[-1.5 0.4]; zl=[-0.16,1.5];  
    case 9
        subdir = '0r-9-700-1-1050-8-700-3-1050';
        fn = 'batchScanned20200620220034.txt';
        xl=[-0.2,1.5]; yl=[-1.5 0.4]; zl=[-0.16,1.5];   
    case 10
        subdir = '0r-9-700-1-1400-8-700-3-1400';
        fn = 'batchScanned20200620220524.txt';
        xl=[-0.2,1.5]; yl=[-1.5 0.4]; zl=[-0.16,1.5]; 
    case 11
        subdir = '0r-9-700-3-1050-8-700-10-1050';
        fn = 'batchScanned20200621232119.txt';
        xl=[-0.2,1.5]; yl=[-1.5 0.4]; zl=[-0.16,1.5]; 
    case 12
        subdir = '0r-9-700-3-1400-8-700-10-1400';
        fn = 'batchScanned20200621231636.txt';
        xl=[-0.2,1.5]; yl=[-1.5 0.4]; zl=[-0.16,1.5]; 
    case 13
        subdir = '1r-9-350-1-700-8-350-3-700';
        fn = 'batchScanned20200620223207.txt';
        xl=[-1.5,0.2]; yl=[-0.4 1.5]; zl=[-0.16,1.5];   
    case 14
        subdir = '1r-9-350-1-1050-8-350-3-1050';
        fn = 'batchScanned20200620222711.txt';
        xl=[-1.5,0.2]; yl=[-0.4 1.5]; zl=[-0.16,1.5];     
    case 15
        subdir = '1r-9-350-3-700-8-350-10-700';
        fn = 'batchScanned20200621215243.txt';
        xl=[-1.5,0.2]; yl=[-0.4 1.5]; zl=[-0.16,1.5];    
    case 16
        subdir = '1r-9-350-3-1050-8-350-10-1050';
        fn = 'batchScanned20200621221312.txt';
        xl=[-1.5,0.2]; yl=[-0.4 1.5]; zl=[-0.16,1.5];   
    case 17
        subdir = '1r-9-520-1-1050-8-520-3-1050';
        fn = 'batchScanned20200620222226.txt';
        xl=[-1.5,0.2]; yl=[-0.4 1.5]; zl=[-0.16,1.5];   
    case 18
        subdir = '1r-9-520-1-1400';
        fn = 'batchScanned20200621214519.txt';
        xl=[-1.8,0.2]; yl=[-0.4 1.5]; zl=[-0.16,1.5]; 
    case 19
        subdir = '1r-9-520-3-1050-8-520-10-1050';
        fn = 'batchScanned20200621223719.txt';
        xl=[-1.5,0.2]; yl=[-0.4 1.5]; zl=[-0.16,1.5];  
    case 20
        subdir = '1r-9-520-3-1400-8-520-10-1400';
        fn = 'batchScanned20200621223045.txt';
        xl=[-1.6,0.2]; yl=[-0.4 1.6]; zl=[-0.16,1.5]; 
    case 21
        subdir = '1r-9-700-1-1050-8-700-3-1050';
        fn = 'batchScanned20200620223726.txt';
        xl=[-1.6,0.2]; yl=[-0.4 1.6]; zl=[-0.16,1.5]; 
    case 22
        subdir = '1r-9-700-1-1400-8-700-3-1400';
        fn = 'batchScanned20200620224153.txt';
        xl=[-1.6,0.2]; yl=[-0.4 1.6]; zl=[-0.16,1.5]; 
    case 23
        subdir = '1r-9-700-3-1050-8-700-10-1050';
        fn = 'batchScanned20200621221916.txt';
        xl=[-1.6,0.2]; yl=[-0.4 1.6]; zl=[-0.16,1.5]; 
    case 24
        subdir = '1r-9-700-3-1400-8-700-10-1400';
        fn = 'batchScanned20200621222412.txt';
        xl=[-1.6,0.2]; yl=[-0.4 1.6]; zl=[-0.16,1.5]; 
    case 25
        subdir = '1r-1-2000-13-3000';
        fn = 'batchScanned20201024142905.txt';
        xl=[-4.5,0.2]; yl=[-0.5 0.5]; zl=[-0.16,1.5]; 
    case 26
        subdir = '1r-1-2000-13-4000';
        fn = 'batchScanned20201024143559.txt';
        xl=[-4.5,0.2]; yl=[-0.5 0.5]; zl=[-0.16,1.5]; 
    case 27
        subdir = '1r-3-2000-13-3000';
        fn = 'batchScanned20201024144943.txt';
        xl=[-4.5,0.2]; yl=[-0.5 0.5]; zl=[-0.16,1.5]; 
    case 28
        subdir = '1r-3-2000-13-4000';
        fn = 'batchScanned20201024144256.txt';
        xl=[-4.5,0.2]; yl=[-0.5 0.5]; zl=[-0.16,1.5]; 
    case 29
        subdir = '1r-13-2000-1-3000';
        fn = 'batchScanned20201024141019.txt';
        xl=[-4.5,0.2]; yl=[-0.5 0.5]; zl=[-0.16,1.5]; 
    case 30
        subdir = '1r-13-2000-1-4000';
        fn = 'batchScanned20201024142110.txt';
        xl=[-4.5,0.2]; yl=[-0.5 0.5]; zl=[-0.16,1.5]; 
    case 31
        subdir = '1r-R-2000-B-3000';
        fn = 'batchScanned20201024153218.txt';
        xl=[-4.5,0.2]; yl=[-0.5 0.5]; zl=[-0.16,1.5]; 
    case 32
        subdir = '1r-R-2000-B-4000';
        fn = 'batchScanned20201024152549.txt';
        xl=[-4.5,0.2]; yl=[-0.5 0.5]; zl=[-0.16,1.5]; 
    case 33
        subdir = '1r-Y-2000-B-3000';
        fn = 'batchScanned20201024151138.txt';
        xl=[-4.5,0.2]; yl=[-0.5 0.5]; zl=[-0.16,1.5]; 
    case 34
        subdir = '1r-Y-2000-B-4000';
        fn = 'batchScanned20201024151836.txt';
        xl=[-4.5,0.2]; yl=[-0.5 0.5]; zl=[-0.16,1.5]; 
    case 35
        subdir = '1r-Y-2000-R-3000';
        fn = 'batchScanned20201024145754.txt';
        xl=[-4.5,0.2]; yl=[-0.5 0.5]; zl=[-0.16,1.5]; 
    case 36
        subdir = '1r-Y-2000-R-4000';
        fn = 'batchScanned20201024150434.txt';
        xl=[-4.5,0.2]; yl=[-0.5 0.5]; zl=[-0.16,1.5]; 
    otherwise % put new data here
        dataname = 'batchScanned20200703172752';
end

% read data
% x1 = -2.5; x2 = 1;
% y1 = -1; y2 = 1.8;
% z1 = -0.3; z2 = 0.9;

[range, angleV, angleH, timestamp] = read_scandata([dir,'\',subdir,'\',fn]); 
range = replace_outlier(range,0.1, 10);

%% show point cloud
if subdir(1)=='0'
    R=0.0387; Dphi=0.0000; Dpsi=0.4000; Dtheta=-1.84;
else
    R=0.1919; Dphi=0.0000; Dpsi=0.4000; Dtheta=-1.74;
end
ps = range2pointsPrecise(range, angleV, angleH, R,Dphi,Dpsi,Dtheta);
figure(1); 
scatter3(ps(:,1),ps(:,2),ps(:,3),1,'.');
az = 30; el = 60; view(az,el)
% xlim([-2,0.3]); ylim([-0.4 2]); 
% zlim([-0.16,0.2])
xlim(xl); ylim(yl); zlim(zl)
xlabel('x'); ylabel('y'); zlabel('z'); 

