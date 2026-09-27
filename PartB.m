%% Aerodinamica, Mecanica de Vol i Orbital - Project: Part B
% Pressure-velocity coupling verification (Projection Method)
% Reference: Topic 2, Slides 61 & 63.

clear; clc;

%% 1. Mesh definition
N = 3;
M = 3;
L = 1.0;
H = 1.0;

mesh = create_mesh(N, M, L, H);

%% 2. Arbitrary predictor velocity field (Slide 63)
% Zero everywhere except a single perturbation: u^p(3,3) = 1.0
u_p = zeros(N + 2, M + 2);
v_p = zeros(N + 2, M + 2);
u_p(3, 3) = 1.0;

u_p = halo_update(u_p);
v_p = halo_update(v_p);

%% 3. Step 1: Integrated divergence of u^p (Slide 61, arrow 1)
d_p = diverg(u_p, v_p, mesh);
print_field(d_p, 'Divergence of u^p (before correction)');

% Check global mass balance (sum of Poisson source terms must be zero)
sum_div_p = sum(d_p(2:N+1, 2:M+1), 'all');
fprintf('Sum of divergence of u^p: %e\n\n', sum_div_p);

%% 4. Step 2: Laplacian matrix (Slide 61, arrow 2)
A = laplacianMatrix(mesh);

%% 5. Step 3: Field to algebraic vector (Slide 61, arrow 3)
b = field2vector(d_p, mesh);

%% 6. Step 4: Solve Poisson linear system: A * p = b (Slide 61, step 4)
p = A \ b;

%% 7. Step 5: Vector to 2D field and halo update (Slide 61, arrow 5)
P = vector2field(p, mesh);

%% 8. Step 6: Gradient of pseudo-pressure & velocity correction (Slide 61, arrow 6)
[gx, gy] = gradient(mesh, P);

u_next = halo_update(u_p - gx);
v_next = halo_update(v_p - gy);

%% 9. Step 7: Divergence of corrected velocity u^{n+1} (Slide 61, arrow 7 & Slide 63)
d_next = diverg(u_next, v_next, mesh);
print_field(d_next, 'Divergence of u^{n+1} (after correction)');

%% 10. Numerical verification summary
max_div_before = max(abs(d_p(2:N+1, 2:M+1)), [], 'all');
max_div_after  = max(abs(d_next(2:N+1, 2:M+1)), [], 'all');

fprintf('Maximum divergence before correction: %e\n', max_div_before);
fprintf('Maximum divergence after correction : %e\n', max_div_after);
