function [xp,yp,xu,yu,xv,yv] = create_mesh(N,M,L,H)
% [xp,yp,xu,yu,xv,yv] = create_mesh(N,M,L,H)
% Builds the uniform staggered (Arakawa-C) mesh coordinate vectors for
% the periodic domain [0,L]x[0,H], including 1 halo layer on each side
% (xp,xu,xv have N+2 entries; yp,yu,yv have M+2 entries; indices 2:end-1
% are the real domain).
% Pressure nodes are defined first; u and v are derived from them by
% shifting half a cell (hx/2 or hy/2) in x and y respectively.
% Input:
%   N,M: number of control volumes in x and y
%   L,H: domain size in x and y
% Output:
%   xp,yp: coordinates of the pressure (cell-centered) nodes
%   xu,yu: coordinates of the u nodes (staggered in x)
%   xv,yv: coordinates of the v nodes (staggered in y)
% Author: Pablo Urioste // Marc Antich // Martí Esquerda

hx = L/N;
hy = H/M;
i = (1:N+2)';
j = (1:M+2)';

xp = (i-1.5)*hx;
yp = (j-1.5)*hy;

xu = xp + hx/2;
yu = yp;

xv = xp;
yv = yp + hy/2;

end
