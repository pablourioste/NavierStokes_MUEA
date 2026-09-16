function [xp,yp,xu,yu,xv,yv] = create_mesh(N,L)
% [xp,yp,xu,yu,xv,yv] = create_mesh(N,L)
% Builds the uniform staggered (Arakawa-C) mesh coordinate vectors for
% the periodic NxN domain [0,L]x[0,L], including 1 halo layer on each
% side (vectors have N+2 entries, indices 2:N+1 are the real domain).
% Pressure nodes are defined first; u and v are derived from them by
% shifting half a cell (h/2) in x and y respectively.
% Input:
%   N: number of control volumes per axis
%   L: domain size (same in x and y)
% Output:
%   xp,yp: coordinates of the pressure (cell-centered) nodes
%   xu,yu: coordinates of the u nodes (staggered in x)
%   xv,yv: coordinates of the v nodes (staggered in y)
% Author: Pablo Urioste // Marc Antich // Martí Esquerda

h = L/N;
i = (1:N+2)';

xp = (i-1.5)*h;
yp = (i-1.5)*h;

xu = xp + h/2;
yu = yp;

xv = xp;
yv = yp + h/2;

end
