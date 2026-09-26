%% Part B: pressure-velocity coupling 
clear
N = 3; M = 3; L = 1; H = 1;

% 0) Arbitrary predictor velocity 
u = zeros(N+2, M+2);
v = zeros(N+2, M+2);
u(3,3) = 1;
u = halo_update(u);
v = halo_update(v);

% 1) Divergence of up (integral over each pressure CV)
d = diverg(u, v, N, M, L, H);
% 2) Laplacian matrix
A = laplacianMatrix(N, M, L, H); % Es podria dir que aquesta funció guarda com els pesos de cada part
% 3) Field -> vector
b = field2vector(d, N, M);
% (solve the system)
p = A\b;
% 4) Vector -> field, and halo (grad_p needs the neighbours)
P = vector2field(p, N, M);
P = halo_update(P);

% 5)-6) Gradient of the pseudo-pressure at the staggered nodes
[gx, gy] = grad_p(P, N, M, L, H);
% Correction: u^{n+1} = up - grad(p)
u_new = u - gx; % Apliquem Ta hemholtz aixi si o si vel diverg =0
v_new = v - gy; % Apliquem Ta hemholtz aixi si o si vel diverg =0
u_new = halo_update(u_new);
v_new = halo_update(v_new);

% 7) Check: divergence of u^{n+1} must be ~0 at ALL nodes (s'ha de cumplir div nula)
d_new = diverg(u_new, v_new, N, M , L, H);
print_field(d)          % before: NOT zero
print_field(d_new)      % after: ~1e-16 everywhere
max(max(abs(d_new(2:N+1, 2:M+1))))
