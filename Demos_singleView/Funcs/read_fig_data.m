function [X,Y,Z]=read_fig_data(filename) 
% Function: Read X,Y,Z data from fig file.
% Input:
%     filename - fig file name. 
% Output:
%     X,Y,Z - coordinates values. 
%    
% Demo:
% [X,Y,Z]=read_fig_data('comparison-112.fig');
% 
% Writen by LI, Shuqing (jay_issaac@163.com), 20191206
%
    h1 = openfig(filename,'invisible'); % open figure
    D1=get(h1,'Children'); %get the handle of the line object
    X=get(D1,'XData'); %get the x data
    Y=get(D1,'YData'); %get the y data
    Z=get(D1,'ZData'); %get the y data
end