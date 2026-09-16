function [u,v] = set_velocity_field(xu,yu,xv,yv,fu,fv)
% [u,v] = set_velocity_field(xu,yu,xv,yv,fu,fv)
% Evaluates the analytic functions fu(x,y), fv(x,y) at the u and v nodes.
% Halo included (halo coordinates are the periodic images).
% Author: Pablo Urioste // Marc Antich // Martí Esquerda

N = length(xu);
M = length(yu);
u = zeros(N,M);
v = zeros(N,M);

for i = 1:N
    for j = 1:M
        u(i,j) = fu(xu(i),yu(j));
        v(i,j) = fv(xv(i),yv(j));
    end
end
end
