function [gx, gy] = gradient(mesh, s)
%% GRADIENT Computes the gradient of a cell-centered scalar field at staggered nodes.
% Author: Pablo Urioste // Marc Antich // Martí Esquerda // Iván Aguilar
%
% Evaluates point derivatives at the staggered velocity locations:
%   gx(i, j) = dp/dx at u-nodes (vertical cell faces)
%   gy(i, j) = dp/dy at v-nodes (horizontal cell faces)

N = mesh.N;
M = mesh.M;

gx = zeros(N+2, M+2);
gy = zeros(N+2, M+2);

for j = 2:M+1
    dy = mesh.yp(j+1) - mesh.yp(j);
    for i = 2:N+1
        dx = mesh.xp(i+1) - mesh.xp(i);
        gx(i, j) = (s(i+1, j) - s(i, j)) / dx;
        gy(i, j) = (s(i, j+1) - s(i, j)) / dy;
    end
end

gx = halo_update(gx);
gy = halo_update(gy);
end