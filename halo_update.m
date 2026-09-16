function U = halo_update(u)
% U = halo_update(u)
% Fills the halo (first and last row/column) of a periodic field by
% copying the values from the opposite side of the real domain. The
% inner nodes are left untouched. Works for non-square fields.
% Input:
%   u: field with a 1-node halo on each side
% Output:
%   U: same field with the halo updated
% Author: Pablo Urioste // Marc Antich // Martí Esquerda

N = size(u,1);
M = size(u,2);
U = zeros(N,M);

% Inner nodes are preserved
for i = 2:N-1
    for j = 2:M-1
        U(i,j) = u(i,j);
    end
end

% Bottom and top halo rows (periodic in y)
for i = 2:N-1
    U(i,1) = u(i,M-1);
    U(i,M) = u(i,2);
end

% Left and right halo columns (periodic in x)
for j = 2:M-1
    U(1,j) = u(N-1,j);
    U(N,j) = u(2,j);
end
end
