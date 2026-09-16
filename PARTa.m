%% Aerodinamica, Mecanica de Vol i Orbital - Project: Part A
% Author: Pablo Urioste // Marc Antich // Martí Esquerda
%Soc noob en github estic provant
% Last version: 16/09/2026 --18:07

main();

function main()
% Entry point: builds a test field and prints it
F = build_test_field(5,5);
print_field(F);
% Perform a halo update on the test field and print the updated field
U = halo_update(F);
print_field(U);

end

function F = build_test_field(N,M)
% Creates an NxM matrix where F(i,j) encodes its own (i,j) position
F = zeros(N,M);
for i = 1:N
    for j = 1:M
        F(i,j) = i + j/10;
    end
end
end

function print_field(u)
% Prints matrix u as a field: x (i) left-to-right, y (j) bottom-to-top
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

function U=halo_update(u)
% Perform a halo update on the matrix u
N=size(u,1);
M=size(u,2);
U=zeros(N,M);

% We preserve the inner part intact
for i=2:N-1
    for j=2:M-1
        U(i,j)=u(i,j);
    end
end

% We modify and make the halo update
for i=2:N-1
    U(i,1)=u(i,M-1); %Bottom
    U(i,M)=u(i,2); %Top
end

for j=2:M-1
    U(1,j)=u(N-1,j); %Left
    U(N,j)=u(2,j); %Right
end
end
