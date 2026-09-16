%% Aerodinamica, Mecanica de Vol i Orbital - Project: Part B
% Author: Pablo Urioste // Marc Antich // Martí Esquerda
% Implement and verify pressure-velocity coupling

main();

function main()
% Entry point: builds the mesh and verifies the pressure-velocity
% coupling subprograms
N = 8;
M = 8;
L = 1;
H = 1;
[xp,yp,xu,yu,xv,yv] = create_mesh(N,M,L,H);

% TODO: verify pressure-velocity coupling (e.g. Poisson solver /
% projection method) once implemented
end
