function b = field2vector(d, mesh)
%% FIELD2VECTOR Converts a 2D scalar field (with halo) into an algebraic column vector.
% Author: Pablo Urioste // Marc Antich // Martí Esquerda // Iván Aguilar
%
% Orders interior cells row by row (bottom to top, left to right):
% (2,2) -> 1, (3,2) -> 2, ..., (N+1, M+1) -> N*M (slide 62)

N = mesh.N;
M = mesh.M;


b = zeros(N*M, 1);
k = 0;

for j = 2:M+1          
    for i = 2:N+1
        k = k + 1;
        %k= (j - 2)*N + (i - 1); Esta formula funciona correctamente pero
        %con contador es más rapido .
        b(k) = d(i, j);
    end
end
end