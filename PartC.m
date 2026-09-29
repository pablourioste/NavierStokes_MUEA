%% Aerodinamica, Mecanica de Vol i Orbital - Project: Part C
% Transient Navier-Stokes Solver (Adams-Bashforth 2 + Projection Method)
% Reference: Topic 2, Slides 65 & 66.

clear; clc; close all;

%% 1. Mesh definition & physical parameters
N = 16;
M = 16;
L = 1.0;
H = 1.0;

nu = 0.01;
t_final = 10;
f_cfl = 0.1;

mesh = create_mesh(N, M, L, H);

% Precomputation of Laplacian matrix (computed once outside the loop)
A = laplacianMatrix(mesh);

% Initialization of velocity fields (Taylor-Green Vortex)
[u1_sym, u2_sym] = get_analytical_solution();
fu = matlabFunction(u1_sym, 'Vars', [sym('x'), sym('y')]);
fv = matlabFunction(u2_sym, 'Vars', [sym('x'), sym('y')]);

[u, v] = set_velocity_field(mesh.xu, mesh.yu, mesh.xv, mesh.yv, fu, fv);
u = halo_update(u);
v = halo_update(v);

t = 0;
step = 0;

Ru_prev = [];
Rv_prev = [];

%% 2. Temporal Loop (Adams-Bashforth 2)
fprintf('Starting transient simulation up to t = %.3f s...\n', t_final);

while t < t_final
    step = step + 1;
    
    % Adaptive time step based on CFL / von Neumann stability (Slide 66)
    [dt, dt_c, dt_d] = time_stability(mesh, u, v, nu, f_cfl);
    
    % Adjust the final time step to reach t_final exactly
    if t + dt > t_final
        dt = t_final - t;
    end

    % Step 1: Evaluate spatial residual R(u) = -C(u) + nu*D(u)
    cu = convective_u(u, v, mesh);
    cv = convective_v(u, v, mesh);
    du = diffusive_u(u, mesh);
    dv = diffusive_v(v, mesh);
    Ru = -cu + nu * du;
    Rv = -cv + nu * dv;

    % Step 2: Temporal predictor (up, vp)
    if step == 1
        % Bootstrap with forward Euler at the first time step
        up = u + dt * Ru;
        vp = v + dt * Rv;
    else
        % Second-order Adams-Bashforth (AB2)
        up = u + dt * (1.5 * Ru - 0.5 * Ru_prev);
        vp = v + dt * (1.5 * Rv - 0.5 * Rv_prev);
    end

    up = halo_update(up);
    vp = halo_update(vp);

    % Step 3: Poisson projection and velocity correction
    [u, v, P, div_max] = project_velocity(up, vp, A, mesh);

    % Step 4: Update state for the next time step
    Ru_prev = Ru;
    Rv_prev = Rv;
    t = t + dt;

    if mod(step, 25) == 0 || t >= t_final
        fprintf('Step %4d | t = %.4f s | dt = %.2e s | Max Div = %.2e\n', ...
            step, t, dt, div_max);
    end
end

fprintf('\nSimulation completed successfully at t = %.4f s (Total steps: %d).\n', t, step);

%% Helper functions
function [u1, u2] = get_analytical_solution()
syms x y;
u1 = cos(2*pi*x)*sin(2*pi*y);
u2 = -sin(2*pi*x)*cos(2*pi*y);
end