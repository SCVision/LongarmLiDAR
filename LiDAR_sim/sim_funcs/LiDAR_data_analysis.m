function [phi_obj,dens_p,dens_n,sin_beta_p,sin_beta_n,...
    objrng_p,objrng_n,objpos_p,objpos_n] = ...
    LiDAR_data_analysis(phi,range_p,range_n,pos,arm_len,hG,max_rang,d_e)
% Function: simulate adaptive scanning of a long-arm lidar in a 2D map (top view).
% Input:
%     phi - array of long-arm orientaion (Kx1), deg.
%     range_p, range_n - array of range corresponding to phi (Kx1).
%     pos - X, Y coordinates of the LiDAR in the map.
%     arm_len - arm length.
%     hG - width of Gaussian
%     max_rang - range for searching objects
%     d_e - precision for line searching
% Output:
%     phi_obj - array of long-arm orientaion (Nx1), deg.
%     dens_p, dens_n - density corresponding to phi_obj (Nx1).
%     sin_beta_p,sin_beta_n - inclination corresponding to phi_obj (Nx1).
%     objrng_p, objrng_n - range corresponding to phi_obj (Nx1).
%     objpos_p,objpos_n - coordinates of points corresponding to phi_obj (Nx2).
% Demo:
% 
% Writen by LIN, Jingyu (linjy02@hotmail.com), 20210128
%

% paramteres
Kdeg2rad = pi/180;
% hG = e_s; % width for density estimation
% max_rang = 3000; % range for searching objects
% d_e = 0.5; % precision for line searching

phi_quan = 0.5; % resolution of phi
dPhi_dens = 4*max(abs(diff(phi))); % angle range of data for density estimation
di_dens = ceil(dPhi_dens/phi_quan); % angle range of data for density estimation
dPhi_dens = di_dens*phi_quan; 
N = 1 + round((phi(end) - phi(1))/phi_quan); % quantum of phi

% recorders for object detection
phi_obj = phi(1) + (di_dens : N-1-di_dens)'*phi_quan;
dens_p = zeros(length(phi_obj),1);     % density 
dens_n = zeros(length(phi_obj),1); 
sin_beta_p = zeros(length(phi_obj),1); % inclination  
sin_beta_n = zeros(length(phi_obj),1); 
objrng_p = zeros(length(phi_obj),1);   % object range
objrng_n = zeros(length(phi_obj),1); 
objpos_p = zeros(length(phi_obj),2);   % object position
objpos_n = zeros(length(phi_obj),2);

% preparing
[pnts_p,pnts_n] = range2points_topview(pos,phi,range_p,range_n,arm_len);
p_densc = 2; % pointer of phi for data near density
p_dens0 = 1; % pointer of phi for first data for density
p_dens1 = 2; % pointer of phi for last data for density
for i = 2 : length(phi)-1
    if phi(i-1) <= phi_obj(1) && phi(i) > phi_obj(1)
        p_densc = i;
    end
    if phi(i) > phi_obj(1) + dPhi_dens
        p_dens1 = i;
        break;
    end
end

% start scanning
x_pos = pos(1); y_pos = pos(2);
for i_dens = 1 : length(phi_obj)
    %%% density estimation %%%
    phi_dens = phi_obj(i_dens); % angle for density
    if phi(p_dens0+1) < phi_dens - dPhi_dens
        p_dens0 = p_dens0 + 1; % move p_dens0 if it is out of range
    end
    if phi(p_dens1) < phi_dens + dPhi_dens
        p_dens1 = p_dens1 + 1; % move p_dens1 if it is out of range
    end
    ps_p = pnts_p(p_dens0:p_dens1,:); % front data for density estimation
    Xc = x_pos + arm_len*cosd(phi_dens);
    Yc = y_pos - arm_len*sind(phi_dens); % down is positive
    [Xo,Yo,dist,dens] = onMap_dens(Xc,Yc,phi_dens+90,ps_p,hG,max_rang,d_e);
    dens_p(i_dens) = dens;     % 1 record density
    objrng_p(i_dens) = dist;   % 2 record object range
    objpos_p(i_dens,:) = [Xo,Yo];
    ps_n = pnts_n(p_dens0:p_dens1,:); % back data for density estimation
    [Xo,Yo,dist,dens] = onMap_dens(Xc,Yc,phi_dens-90,ps_n,hG,max_rang,d_e);
    dens_n(i_dens) = dens;     % 1 record density
    objrng_n(i_dens) = dist;   % 2 record object range
    objpos_n(i_dens,:) = [Xo,Yo];
    
    %%% object detection & step update %%%
    if phi(p_densc) <= phi_dens
        p_densc = p_densc + 1; % move p_densc if it is out of range
    end
    phi_step_scan = (phi(p_densc)-phi(p_densc-1)); % real phi step
    % front data
    e_m_p = 1/dens_p(i_dens); % point space
    sin_beta_p(i_dens) = objrng_p(i_dens) * phi_step_scan*Kdeg2rad * dens_p(i_dens);
    % back data
    e_m_n = 1/dens_n(i_dens); % point space
    sin_beta_n(i_dens) = objrng_n(i_dens) * phi_step_scan*Kdeg2rad * dens_n(i_dens);
end % for i
