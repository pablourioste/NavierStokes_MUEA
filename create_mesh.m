function mesh = create_mesh(N, M, L, H)
% CREATE_MESH Builds uniform staggered (Arakawa-C) mesh structure with 1 halo layer.
%
% Usage:
%   mesh = create_mesh(N, M, L, H);
%   [xp, yp, xu, yu, xv, yv] = create_mesh(N, M, L, H);
%
% Output:
%   mesh - Structure containing node coordinates, control volume surfaces, and volumes.
%   xp, yp, xu, yu, xv, yv - Node coordinate vectors (if called with multiple outputs).

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

% 2. Coordinate matrices for node locations (oriented i=x, j=y)
[XP, YP] = ndgrid(xp, yp);
[XU, YU] = ndgrid(xu, yu);
[XV, YV] = ndgrid(xv, yv);

% 3. Geometrical properties (assuming unit depth W = 1)
% Convention: first index is x (N+2), second index is y (M+2)
Sx = hy * ones(N+2, M+2); % Surfaces perpendicular to x-axis
Sy = hx * ones(N+2, M+2); % Surfaces perpendicular to y-axis

% Control volumes
Vp = hx * hy * ones(N+2, M+2); % Pressure control volume
Vu = hx * hy * ones(N+2, M+2); % u-velocity control volume
Vv = hx * hy * ones(N+2, M+2); % v-velocity control volume

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