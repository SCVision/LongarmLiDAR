function [range, angleV, angleH, timestamp] = read_scandata(txtfile) 
% Function: import raw data scanned by our 3D Lidar.
%     txtfile - data file name (txt format).
% Output:
%     range - range data (H*V). 
%     angleV - vertical angles theta (V*1).
%     angleH - horizontal angles phi (H*1). 
%     timestamp - timestamp of rows (H*1)
% Demo:
% [range, angleV, angleH, timestamp] = read_scandata('Scanned1.txt'); 
% [V, H]=meshgrid(angleV, angleH);
% figure(1); mesh(V,H,range)
% figure(2); plot(timestamp)
% 
% Writen by LIN, Jingyu (linjy02@hotmail.com), 20200127
%
rawdata = importdata(txtfile);
angleV = rawdata(1, 4:end)';
timestamp = rawdata(2:end, 2);
angleH = rawdata(2:end, 3);
range = rawdata(2:end, 4:end);

