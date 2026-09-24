
%% Aerodinamica, Mecanica de Vol i Orbital - Project: Operator Verification
% Author: Pablo Urioste // Marc Antich // Martí Esquerda // Iván Aguilar
% Verifies spatial convective and diffusive operators using the Method of
% Manufactured Solutions (MMS) with the Taylor-Green vortex benchmark.

% 1. Symbolic velocity field and exact transport operators
syms x y;
u1 = cos(2*pi*x)*sin(2*pi*y);
u2 = -sin(2*pi*x)*cos(2*pi*y);

% Evaluate convective and diffusive terms symbolically:
conv_sym_u = diff(u1^2, x) + diff(u1 * u2, y);
conv_sym_v = diff(u1 * u2, x) + diff(u2^2, y);
diff_sym_u = diff(u1, x, 2) + diff(u1, y, 2);
diff_sym_v = diff(u2, x, 2) + diff(u2, y, 2);


% 2. Vector of grid resolutions (N x N)
N_vec = [10, 20, 40, 80, 160];

% Preallocate maximum error vectors
err_max_cu = zeros(size(N_vec));
err_max_cv = zeros(size(N_vec));
err_max_du = zeros(size(N_vec));
err_max_dv = zeros(size(N_vec));

for k = 1:length(N_vec)
    N = N_vec(k);
    mesh = create_mesh(N, N, 1, 1);
    
    % Node coordinates:
    nodes_u = zeros(mesh.N+2, mesh.M+2, 2);
    nodes_u(:,:,1) = mesh.XU;
    nodes_u(:,:,2) = mesh.YU;

    nodes_v = zeros(mesh.N+2, mesh.M+2, 2);
    nodes_v(:,:,1) = mesh.XV;
    nodes_v(:,:,2) = mesh.YV;

    % Analytical function evaluations:
    ua = evaluation(u1, nodes_u);
    va = evaluation(u2, nodes_v);

    cu_analytic = evaluation(conv_sym_u, nodes_u);
    cv_analytic = evaluation(conv_sym_v, nodes_v);
    du_analytic = evaluation(diff_sym_u, nodes_u);
    dv_analytic = evaluation(diff_sym_v, nodes_v);
    
    % Numerically computed operators:
    cu = convective_u(ua, va, mesh);
    cv = convective_v(ua, va, mesh);
    du = diffusive_u(ua, mesh);
    dv = diffusive_v(va, mesh);

    % Maximum error over interior nodes:
    err_max_cu(k) = max(abs(cu_analytic(2:end-1, 2:end-1) - cu(2:end-1, 2:end-1)), [], 'all');
    err_max_cv(k) = max(abs(cv_analytic(2:end-1, 2:end-1) - cv(2:end-1, 2:end-1)), [], 'all');
    err_max_du(k) = max(abs(du_analytic(2:end-1, 2:end-1) - du(2:end-1, 2:end-1)), [], 'all');
    err_max_dv(k) = max(abs(dv_analytic(2:end-1, 2:end-1) - dv(2:end-1, 2:end-1)), [], 'all');
end

% Characteristic spatial step h = 1/N
h = 1 ./ N_vec;

figure('Color', 'w');

% Plot numerical errors versus h = 1/N
loglog(h, err_max_cu, '-o', 'LineWidth', 2, 'DisplayName', 'Convective U'); hold on;
loglog(h, err_max_cv, '-s', 'LineWidth', 2, 'DisplayName', 'Convective V');
loglog(h, err_max_du, '-^', 'LineWidth', 2, 'DisplayName', 'Diffusive U');
loglog(h, err_max_dv, '-d', 'LineWidth', 2, 'DisplayName', 'Diffusive V');

% Theoretical reference slopes as a function of h
% Under grid refinement, h decreases and error drops proportionally to h^p
ref_1st = err_max_du(1) * (h ./ h(1));       % O(h^1) - 1st order
ref_2nd = err_max_du(1) * (h ./ h(1)).^2;   % O(h^2) - 2nd order

loglog(h, ref_1st, '--k', 'LineWidth', 1.2, 'DisplayName', '\mathcal{O}(h^1) (1^{st} order)');
loglog(h, ref_2nd, '-.k', 'LineWidth', 1.2, 'DisplayName', '\mathcal{O}(h^2) (2^{nd} order)');

hold off;
grid on; box on;
xlabel('Grid spacing h = 1/N', 'FontSize', 12);
ylabel('Maximum Error ||e||_{\infty}', 'FontSize', 12);
title('Spatial Grid Convergence Analysis of Operators', 'FontSize', 14);
legend('Location', 'southeast', 'FontSize', 10);