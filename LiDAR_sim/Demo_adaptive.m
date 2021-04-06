% simulation of long-arm LiDAR in adaptive mode
path(path,'.\sim_funcs')
adp = 1; % switch adaptive scanning

% parameters
arm_len = 100;
samples = 240;
show_step = 3;
d_s = 500; % working distance
e_s = 2*pi/samples*d_s; % expected point space
show_p = 1; show_n = 1; % time for showing£¨ms£©
map = double(imread('scene_adaptive.bmp'));
max_rang = 1500; % max range of lidar
pos = [1500, 1500]; % LiDAR position

% scanning
[phi,range_p,range_n, phi_obj,dens_p,dens_n, ...
    sin_beta_p,sin_beta_n,objrng_p,objrng_n,objpos_p,objpos_n] ...
    = LiDAR_Scan_adaptive(map,pos,arm_len,e_s,d_s,max_rang,...
    adp,show_step,show_p,show_n);

%% scanning results
[pnts_p,pnts_n] = range2points_topview(pos,phi,range_p,range_n,arm_len);
if adp>0
    % adaptive mode
    figure(10);
    plot(pnts_p(:,1),size(map,1)-pnts_p(:,2),'.r',...
        pnts_n(:,1),size(map,1)-pnts_n(:,2),'.b')
    xlim([0,size(map,2)]);ylim([0,size(map,1)]);
    legend('p','n'); title('points - adaptive')
    figure(11); 
    plot(phi_obj,dens_p,'r',phi_obj,dens_n,'b')
    ylim([0,0.12]); xlim([0 360])
    legend('p','n'); title('density - adaptive')
else
    % fixed mode
    figure(20);
    plot(pnts_p(:,1),size(map,1)-pnts_p(:,2),'.r',...
        pnts_n(:,1),size(map,1)-pnts_n(:,2),'.b')
    xlim([0,size(map,2)]);ylim([0,size(map,1)]);
    legend('p','n'); title('points - fixed step')
    figure(21); 
    plot(phi_obj,dens_p,'r',phi_obj,dens_n,'b')
    ylim([0,0.12]); xlim([0 360])
    legend('p','n'); title('density - fixed step')
end

