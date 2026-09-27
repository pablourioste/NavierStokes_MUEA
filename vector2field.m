function d = vector2field(b, mesh)
% VECTOR2FIELD Converts an algebraic column vector into a 2D scalar field (with halo).
% Author: Pablo Urioste // Marc Antich // Martí Esquerda // Iván Aguilar
%
% Inverse of field2vector: places vector elements back into interior cells
% row by row (bottom to top, left to right) and updates periodic halos.

N = mesh.N;
M = mesh.M;

d = zeros(N+2, M+2);
k = 0;

for j = 2:M+1
    for i = 2:N+1
        k = k + 1;
        d(i, j) = b(k);
    end
end

d = halo_update(d);
end