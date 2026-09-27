clear; clc; close all;

%% 1. Domain Parameters & Staggered Mesh Generation
N = 16;
M = 16;
L = 1.0;
H = 1.0;

mesh = create_mesh(N, M, L, H);

%% 2. Mesh Visualization
plot_staggered_mesh(mesh);

%% 3. Basic Building Blocks Verification
N_demo = 4;
M_demo = 4;
F_test = build_test_field(N_demo + 2, M_demo + 2);
print_field(F_test, 'Original Test Field (Before Halo Update)');

F_halo = halo_update(F_test);
print_field(F_halo, 'Field with Updated Halos (Periodic Boundaries)');

%% 4. Staggered Velocity Field Initialization
[u1_sym, u2_sym] = get_analytical_solution();
fu = matlabFunction(u1_sym, 'Vars', [sym('x'), sym('y')]);
fv = matlabFunction(u2_sym, 'Vars', [sym('x'), sym('y')]);

[u, v] = set_velocity_field(mesh.xu, mesh.yu, mesh.xv, mesh.yv, fu, fv);
u = halo_update(u);
v = halo_update(v);

%% 5. Spatial Operator Verification
N_vec = [10, 20, 40, 80, 160;
         100, 200, 400, 800, 1600];
verification_results = verification(u1_sym, u2_sym, N_vec);

% Repeat verification and see what happens when we multiply function u1_sym
% by 4.
u1_scaled= u1_sym*4;
u2_scaled= u2_sym;
scaled_verification_results = verification(u1_scaled, u2_scaled, N_vec);


%% Functions
function [u1, u2] = get_analytical_solution()
syms x y;
u1 = cos(2*pi*x)*sin(2*pi*y);
u2 = -sin(2*pi*x)*cos(2*pi*y);
end

function F = build_test_field(N, M)
F = zeros(N, M);
for i = 2:N-1
    for j = 2:M-1
        F(i, j) = i + j / 10;
    end
end
end
