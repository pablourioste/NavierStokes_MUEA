function print_field(F,name)
% print_field(u)
% Prints a 2D field using the project convention: first index (i) is the
% horizontal position (left to right), second index (j) is the vertical
% position (bottom to top). Halo nodes are printed as well.
% Input:
%   u: field to print (matrix, halo included)
% Author: Pablo Urioste // Marc Antich // Martí Esquerda

if nargin < 2 || isempty(name)
    name = 'Field';
end

[M, N] = size(F);
fprintf('\n================ %s (%dx%d) ================\n', name, M, N);


if (N<1000 && M<1000)
    for j = M:-1:1
        fprintf('j=%2d | ', j);
        for i = 1:N
            fprintf('%+6.2f ', F(i, j));
        end
        fprintf('\n');
    end
else
    fprintf('Matrix too big to print.\n');
end
fprintf('===================================================\n\n');
end



