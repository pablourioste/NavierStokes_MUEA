% function [xp,yp,xu,yu,xv,yv] = create_mesh(N,M,L,H)
% % [xp,yp,xu,yu,xv,yv] = create_mesh(N,M,L,H)
% % Builds the uniform staggered (Arakawa-C) mesh coordinate vectors for
% % the periodic domain [0,L]x[0,H], including 1 halo layer on each side
% % (xp,xu,xv have N+2 entries; yp,yu,yv have M+2 entries; indices 2:end-1
% % are the real domain).
% % Pressure nodes are defined first; u and v are derived from them by
% % shifting half a cell (hx/2 or hy/2) in x and y respectively.
% % Input:
% %   N,M: number of control volumes in x and y
% %   L,H: domain size in x and y
% % Output:
% %   xp,yp: coordinates of the pressure (cell-centered) nodes
% %   xu,yu: coordinates of the u nodes (staggered in x)
% %   xv,yv: coordinates of the v nodes (staggered in y)
% % Author: Pablo Urioste // Marc Antich // Martí Esquerda
% 
% hx = L/N;
% hy = H/M;
% i = (1:N+2)';
% j = (1:M+2)';
% 
% xp = (i-1.5)*hx;
% yp = (j-1.5)*hy;
% 
% xu = xp + hx/2;
% yu = yp;
% 
% xv = xp;
% yv = yp + hy/2;
% 
% end


function mesh = create_mesh(N, M, L, H)
% CREATE_MESH Builds uniform staggered (Arakawa-C) mesh structure with 1 halo layer.
%
% Output:
%   mesh - Structure containing node coordinates, control volume surfaces, and volumes.

hx = L / N;
hy = H / M;

i = (1:N+2)';
j = (1:M+2)';

% 1. Node coordinate vectors
xp = (i - 1.5) * hx;
yp = (j - 1.5) * hy;

xu = xp + hx / 2;
yu = yp;

xv = xp;
yv = yp + hy / 2;

% 2. Meshgrid matrices for node locations
[XP, YP] = meshgrid(xp, yp);
[XU, YU] = meshgrid(xu, yu);
[XV, YV] = meshgrid(xv, yv);

% 3. Geometrical properties (assuming unit depth W = 1)
% Control volume cell face surfaces
Sx = hy * ones(M+2, N+2); % Surfaces perpendicular to x-axis
Sy = hx * ones(M+2, N+2); % Surfaces perpendicular to y-axis

% Control volumes
Vp = hx * hy * ones(M+2, N+2); % Pressure control volume
Vu = hx * hy * ones(M+2, N+2); % u-velocity control volume
Vv = hx * hy * ones(M+2, N+2); % v-velocity control volume

% Store in output struct
mesh.N = N; mesh.M = M;
mesh.L = L; mesh.H = H;
mesh.hx = hx; mesh.hy = hy;

mesh.xp = xp; mesh.yp = yp;
mesh.xu = xu; mesh.yu = yu;
mesh.xv = xv; mesh.yv = yv;

mesh.XP = XP; mesh.YP = YP;
mesh.XU = XU; mesh.YU = YU;
mesh.XV = XV; mesh.YV = YV;

mesh.Sx = Sx;
mesh.Sy = Sy;

mesh.Vp = Vp;
mesh.Vu = Vu;
mesh.Vv = Vv;
end