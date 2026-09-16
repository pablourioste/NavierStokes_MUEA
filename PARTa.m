%% Aerodinamica, Mecanica de Vol i Orbital - Project: Part A
% Author: Pablo Urioste // Marc Antich // Martí Esquerda
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
