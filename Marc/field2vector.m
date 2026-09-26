function b = field2vector(d, N, M)
% b = field2vector(d, N, M)
% Converts a centred scalar field (with halo) into an algebraic vector.
% Only interior nodes are stored, ordered row by row (bottom to top,
% left to right), as in slide 62: (2,2)->1, (3,2)->2, ..., (N+1,M+1)->N*M
%   d    : centred scalar field (N+2)x(M+2)
%   N, M : number of interior cells in x and y
%   b    : column vector (N*M)x1

b = zeros(N*M, 1);
k = 0;                      % contador: posición en el vector

for j = 2:M+1               % filas: de abajo arriba
    for i = 2:N+1           % dentro de cada fila: de izquierda a derecha
        k = k + 1;          % siguiente posición
        b(k) = d(i,j);
    end
end

end
