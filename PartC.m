%% Aerodinamica, Mecanica de Vol i Orbital - Project: Part C
% Author: Pablo Urioste // Marc Antich // Martí Esquerda
% Implement and verify time integration

main();

function main()
% Entry point: builds the mesh and verifies the time integration
% subprograms
N = 8;
M = 8;
L = 1;
H = 1;
[xp,yp,xu,yu,xv,yv] = create_mesh(N,M,L,H);

% TODO: verify time integration scheme once implemented
end
