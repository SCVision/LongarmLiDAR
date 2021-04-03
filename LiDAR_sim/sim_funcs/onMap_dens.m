function [Xo,Yo,dist,dens] = onMap_dens(Xc,Yc,theta,pnts,hG,max_rang,d_e)
% Function: search object at a beam of laser.
% Input:
%     Xc,Yc - starting point of laser in the map.
%     theta - direction of laser, deg.
%     pnts - point set for density estimation (Nx2).
%     hG - width of Gaussian
%     max_rang - range for searching objects
%     d_e - precision of range.
% Output:
%     Xo,Yo - position of the object being found.
%     dist - distance of the object being found.
%     dens - density of points on the object.
% 
% Writen by LIN, Jingyu (linjy02@hotmail.com), 20210121
%

if nargin < 7, d_e = 0.2; end  % case: val = func()
% rng = sqrt(pnts*pnts'); % points to ranges
% max_rang = max(rng(:)); 
dens = 0;
dist = max_rang;
Xo = Xc;
Yo = Yc;

% start searching
for t = 1:d_e:max_rang
    X = Xc + t*cosd(theta);
    Y = Yc - t*sind(theta); % down is positive
    pnt_dens = KDE_2D([X,Y], pnts, hG);
    if  dens < pnt_dens
        dens = pnt_dens;
        dist = t;
        Xo = X;
        Yo = Y;
    end
end
