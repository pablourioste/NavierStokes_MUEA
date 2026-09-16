function print_field(u)
% print_field(u)
% Prints a 2D field using the project convention: first index (i) is the
% horizontal position (left to right), second index (j) is the vertical
% position (bottom to top). Halo nodes are printed as well.
% Input:
%   u: field to print (matrix, halo included)
% Author: Pablo Urioste // Marc Antich // Martí Esquerda

N = size(u,1);
M = size(u,2);
fprintf('Matrix u has %d rows and %d columns.\n', N, M);

if (N<100 && M<100)
    for j = M:-1:1
        fprintf('j=%d ', j);
        for i = 1:N
            fprintf('%+8.4f ', u(i,j));
        end
        fprintf('\n');
    end
else
    fprintf('Matrix too big to print.\n');
end
end
