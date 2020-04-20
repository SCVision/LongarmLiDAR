path(path,'..\Funcs')

%% read data
[range, angleV, angleH, timestamp] = ...
    read_scandata('data\batchScanned20191127093702.txt'); 
%% remove outlier
% [range_processed] = remove_min_outlier(range, 0.05);
range_processed = medfilt2(range,[5 5]);
R = 0.1870473368; Dtheta = -2.3113893118;
% R = 0.2009314818; Dtheta = -2.6846063503; % mismatching
ps = range2points(range_processed, angleV, angleH, R, Dtheta);

figure(20); 
scatter3(ps(:,1),ps(:,2),ps(:,3),1,'.'); 
xlabel('x'); ylabel('y'); zlabel('z'); 
xlim([-1.5,2.1]); ylim([-2 5]); zlim([-0.5,3.5])

%% crop and show
x1 = -1; x2 = 1.8;
y1 = -1.5; y2 = 2;
z1 = -0.5; z2 = 2;
box = [x1, y1, z1; x2, y2, z2];
[psin, psout] = crop_pointcloud(ps, box);

fig3d = figure(21);
pnts = scatter3(psin(:,1),psin(:,2),psin(:,3),1,'.');
xlabel('x'); ylabel('y'); zlabel('z');
xlim([-1,2]); ylim([-1.5 2]); zlim([-0.5,2])

%% spinning only
for i = 0:360
    az = i; el = 40;
    %view(az,el)
    camorbit(1,0);
    pause(0.02)
end

%% save spinning video
% fn = 'pointcloud_demo'; % video file name
% profile = 'MPEG-4';     % video profile
% v = VideoWriter(fn,profile); % video object
% v.FrameRate = 30;
% v.Quality = 85;
% open(v)
% for i = 0:360
%     az = i; el = 40;
%     view(az,el)
%     frame = getframe(fig3d);
%     writeVideo(v,frame);
% end
% close(v)