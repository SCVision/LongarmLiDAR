function [phi,range_p,range_n, phi_obj,dens_p,dens_n, ...
    sin_beta_p,sin_beta_n,objrng_p,objrng_n,objpos_p,objpos_n] ...
    = LiDAR_Scan_adaptive(map,pos,arm_len,e_s,d_s,max_rang, ...
    adp,show_step,show_p,show_n)
% Function: simulate adaptive scanning of a long-arm lidar in a 2D map (top view).
% Input:
%     map - array of 2D map (WxH), right & down are positive. 
%           The map defines distance unit.
%     pos - X, Y coordinates of the LiDAR in the map.
%     arm_len - arm length.
%     e_s - expected point space.
%     d_s - working distance.
%     max_rang - range for searching objects
%     adp - switch adaprive scanning
%     show_step - laser is shown erery this step
%     show_p, show_n - time interval when showing laser (ms). 
%                  0 means no showing.
% Output:
%     phi - array of long-arm orientaion (Kx1), deg.
%     range_p, range_n - array of range corresponding to phi (Kx1).
%     phi_obj - array of long-arm orientaion (Nx1), deg.
%     dens_p, dens_n - density corresponding to phi_obj (Nx1).
%     sin_beta_p,sin_beta_n - inclination corresponding to phi_obj (Nx1).
%     objrng_p, objrng_n - range corresponding to phi_obj (Nx1).
%     objpos_p,objpos_n - coordinates of points corresponding to phi_obj (Nx2).
% Demo:
% 
% Writen by LIN, Jingyu (linjy02@hotmail.com), 20210111
%

% paramteres
Kdeg2rad = pi/180;
% max_rang = 3000; % range for searching objects
d_e = 0.5; % precision for line searching
sin_beta_min = 0.3; % no object threshold: 20 deg
phi_s = e_s/d_s/Kdeg2rad; % initial step of phi
phi_s_curr = phi_s; % current step of phi
hG = max_rang*phi_s*Kdeg2rad; % width of Gaussian
% hG = e_s; % width for density estimation
phi_z = 0.2*phi_s; % resolution of phi
dPhi_dens = 6*phi_s; % angle range of data for density estimation
di_dens = ceil(dPhi_dens/phi_z); % angle range of data for density estimation
dPhi_dens = di_dens*phi_z; 
phi_start = -dPhi_dens; % starting phi
phi_end = 360 + dPhi_dens;% ending phi
N = 1 + round((phi_end - phi_start)/phi_z); % quantum of phi

% recorders for scanning
phi_data = zeros(N,1); % long-arm orientaions, deg.
range_p = zeros(N,1); range_n = zeros(N,1); % raw range data
pnts_p = zeros(N,2); pnts_n = zeros(N,2);   % points from range data
p_densc = 2; % pointer of phi for data near density
p_dens0 = 1; % pointer of phi for first data for density
p_phi = 0;   % pointer of phi for scanning

% recorders for object detection
phi_obj = phi_start + (0:N-1)'*phi_z;
dens_p = zeros(N,1); dens_n = zeros(N,1); % density 
sin_beta_p = zeros(N,1); sin_beta_n = zeros(N,1); % inclination  
objrng_p = zeros(N,1); objrng_n = zeros(N,1); % object range
objpos_p = zeros(N,2); objpos_n = zeros(N,2); % object position
phi_step = ones(N,1) * phi_s; % scanning step

% start scanning
x_pos = pos(1); y_pos = pos(2);
phi_lastscan = phi_start - 360;
if show_p>0 || show_n>0
    figure(100); imagesc(map); colormap(gray(256))
    hold on
end
for i_phi = 1:N
%     phi_curr = phi_start + (i_phi-1)*phi_quan;
    phi_curr = phi_obj(i_phi);
    % 1. scanning
    if phi_curr - phi_lastscan >= phi_s_curr
        p_phi = p_phi + 1;
        phi_data(p_phi) = phi_curr; % 1.1 record angle 
        phi_lastscan = phi_curr;
        phi_s_curr = phi_step(i_phi); % update scanning step
        x_phi = x_pos + arm_len*cosd(phi_curr);
        y_phi = y_pos - arm_len*sind(phi_curr); % down is positive
        % front scannning
        [x_p,y_p,dist] = onMap_laser(x_phi,y_phi,phi_curr+90,map,d_e);
        range_p(p_phi) = dist;     % 1.2 record range
        pnts_p(p_phi,:) = [x_p,y_p];
        if show_p>0 && rem(p_phi,show_step)==1
            figure(100); plot([x_phi x_p],[y_phi y_p],'r-')
            pause(show_p/1000)
            drawnow
        end
        % back scanning
        [x_n,y_n,dist] = onMap_laser(x_phi,y_phi,phi_curr-90,map,d_e);
        range_n(p_phi) = dist;     % 1.2 record range
        pnts_n(p_phi,:) = [x_n,y_n];
        if show_n>0 && rem(p_phi,show_step)==1
            figure(100); plot([x_phi x_n],[y_phi y_n],'b-')
            pause(show_n/1000)
            drawnow
        end
%         % convert range to point
%         [ps_p,ps_n] = range2points_topview(pos,phi_curr,...
%             range_p(p_phi),range_n(p_phi),arm_len);
%         pnts_p(p_phi,:) = ps_p;
%         pnts_n(p_phi,:) = ps_n;
    end
    
    if i_phi > di_dens*2
        % 2. density estimation
        i_dens = i_phi - di_dens; % angle index for density
%         phi_dens = phi_start + i_dens*phi_quan; % angle for density
%         phi_dens = phi_curr - dPhi_dens; % angle for density
        phi_dens = phi_obj(i_dens); % angle for density
        if phi_data(p_dens0) < phi_dens - dPhi_dens
            p_dens0 = p_dens0 + 1; % move p_dens0 if it is out of range
        end
        ps_p = pnts_p(p_dens0:p_phi,:); % front data for density estimation
        Xc = x_pos + arm_len*cosd(phi_dens);
        Yc = y_pos - arm_len*sind(phi_dens); % down is positive
        [Xo,Yo,dist,dens] = onMap_dens(Xc,Yc,phi_dens+90,ps_p,hG,max_rang,d_e);
        dens_p(i_dens) = dens;     % 2.1 record density
        objrng_p(i_dens) = dist;   % 2.2 record object range
        objpos_p(i_dens,:) = [Xo,Yo];
        ps_n = pnts_n(p_dens0:p_phi,:); % back data for density estimation
        [Xo,Yo,dist,dens] = onMap_dens(Xc,Yc,phi_dens-90,ps_n,hG,max_rang,d_e);
        dens_n(i_dens) = dens;     % 2.1 record density
        objrng_n(i_dens) = dist;   % 2.2 record object range
        objpos_n(i_dens,:) = [Xo,Yo];
        
        % 3. object detection & step update
%         phi_step_scan = (phi_data(p_phi)-phi_data(p_dens0))...
%             / (p_phi-p_dens0); % real phi step
        while phi_data(p_densc) <= phi_dens
            p_densc = p_densc + 1; % move p_densc if it is out of range
        end
        % real phi step       
        phi_step_scan = (phi_data(p_densc)-phi_data(p_densc-1)); 
        % front data
        e_m_p = 1/dens_p(i_dens); % point space
        sin_beta_p(i_dens) = objrng_p(i_dens) * phi_step_scan*Kdeg2rad * dens_p(i_dens);
        alpha_p = atand(objrng_p(i_dens)/arm_len);
        i_phi_q = i_dens + round(2*alpha_p/phi_z);
        if adp>0 && e_m_p>e_s && i_phi_q<=N ...
                && sin_beta_p(i_dens)>sin_beta_min && sin_beta_p(i_dens)<1
            % object detected
            beta_p = asind(sin_beta_p(i_dens));
            phi_step_new = e_s * sind(beta_p+180-2*alpha_p)...
                / objrng_p(i_dens) / Kdeg2rad;
            if phi_step(i_phi_q) > phi_step_new
                phi_step(i_phi_q) = max(phi_step_new,phi_z);
            end
        end
        % back data
        e_m_n = 1/dens_n(i_dens); % point space
        sin_beta_n(i_dens) = objrng_n(i_dens) * phi_step_scan*Kdeg2rad * dens_n(i_dens);
        alpha_n = atand(objrng_n(i_dens)/arm_len);
        i_phi_q = i_dens + round((360-2*alpha_n)/phi_z);
        if adp>0 && e_m_n>e_s && i_phi_q<=N ...
                && sin_beta_n(i_dens)>sin_beta_min && sin_beta_n(i_dens)<1
            % object detected
            beta_n = asind(sin_beta_n(i_dens));
            phi_step_new = e_s * sind(beta_n+180-2*alpha_n)...
                / objrng_n(i_dens) / Kdeg2rad;
            if phi_step(i_phi_q) > phi_step_new
                phi_step(i_phi_q) = max(phi_step_new,phi_z);
            end
        end
    end % if i_phi is valid for density estimation
    
    % 4. revolve
    if phi_s_curr > phi_step(i_phi)
        phi_s_curr = phi_step(i_phi);
    end
end % for i
if show_p>0 || show_n>0
    figure(100); hold off
end
phi = phi_data(1:p_phi);
range_p = range_p(1:p_phi);
range_n = range_n(1:p_phi); 
phi_obj = phi_obj(1+di_dens:N-di_dens);
dens_p = dens_p(1+di_dens:N-di_dens);
dens_n = dens_n(1+di_dens:N-di_dens);
sin_beta_p = sin_beta_p(1+di_dens:N-di_dens);
sin_beta_n = sin_beta_n(1+di_dens:N-di_dens);
objrng_p = objrng_p(1+di_dens:N-di_dens);
objrng_n = objrng_n(1+di_dens:N-di_dens);
objpos_p = objpos_p(1+di_dens:N-di_dens);
objpos_n = objpos_n(1+di_dens:N-di_dens);