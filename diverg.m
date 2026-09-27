function d = diverg(u, v, mesh)
% DIVERG Evaluates the integrated velocity divergence (net volume flux) on pressure CVs.
% Author: Pablo Urioste // Marc Antich // Martí Esquerda // Iván Aguilar
%
% Returns the net flux out of each control volume:
% This directly matches the RHS vector b for the Poisson solver: A*p = b.

N = mesh.N;
M = mesh.M;

d = zeros(N+2, M+2);

for j = 2:M+1
    for i = 2:N+1
        % Cell face dimensions (surface areas)
        dx = mesh.xu(i) - mesh.xu(i-1);
        dy = mesh.yv(j) - mesh.yv(j-1);

        % Staggered face velocities
        ue = u(i, j);
        uw = u(i-1, j);
        vn = v(i, j);
        vs = v(i, j-1);

        % Net volume flux (integrated divergence across faces)
        d(i, j) = (ue - uw) * dy + (vn - vs) * dx;
    end
end
end