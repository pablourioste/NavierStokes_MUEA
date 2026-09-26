function d = vector2field(b, N, M)
% d = vector2field(b, N, M)
% Converts an algebraic vector back into a centred scalar field.
% Inverse of field2vector: same ordering (slide 62).
% Halo is left at 0: call halo_update afterwards.
%   b    : column vector (N*M)x1
%   N, M : number of interior cells in x and y
%   d    : centred scalar field (N+2)x(M+2)

d = zeros(N+2, M+2);                    % matrix with halo, full of zeros
k = 0;

for j = 2:M+1               % same order as field2vector!
    for i = 2:N+1
        k = k + 1;
        d(i,j) = b(k);                % the assignment, reversed
    end
end

end
