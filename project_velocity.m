function [u_next, v_next, P, max_div] = project_velocity(u_p, v_p, A, mesh)
%% PROJECT_VELOCITY Solves the Poisson equation and projects velocity to divergence-free space.
% Author: Pablo Urioste // Marc Antich // Martí Esquerda // Iván Aguilar
%
% Inputs:
%   u_p   : predictor u-velocity field (size (N+2)x(M+2), halos updated)
%   v_p   : predictor v-velocity field (size (N+2)x(M+2), halos updated)
%   A     : precomputed discrete Laplacian matrix (from laplacianMatrix)
%           or a pre-factorized decomposition object: dA = decomposition(A)
%   mesh  : mesh structure from create_mesh
%
% Outputs:
%   u_next  : corrected divergence-free u-velocity (halos updated)
%   v_next  : corrected divergence-free v-velocity (halos updated)
%   P       : pseudo-pressure field (with halos updated)
%   max_div : (optional) maximum residual divergence in interior cells

% 1. Step 1: Integrated divergence of the predictor field u^p (RHS source term)
d_p = diverg(u_p, v_p, mesh);

% 2. Step 2: Convert 2D divergence field into algebraic column vector b
b = field2vector(d_p, mesh);

% 3. Step 3: Solve Poisson linear system: A * p = b
p = A \ b;

% 4. Step 4: Convert 1D solution back to 2D scalar field (vector2field updates halos)
P = vector2field(p, mesh);

% 5. Step 5: Gradient of pseudo-pressure at staggered cell faces
[gx, gy] = gradient(mesh, P);

% 6. Step 6: Velocity correction (Helmholtz decomposition) and halo update
u_next = halo_update(u_p - gx);
v_next = halo_update(v_p - gy);

% 7. Step 7: (Optional diagnostic) Check residual divergence
if nargout > 3
    N = mesh.N;
    M = mesh.M;
    d_next = diverg(u_next, v_next, mesh);
    max_div = max(abs(d_next(2:N+1, 2:M+1)), [], 'all');
end

end