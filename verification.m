function [results, h_fig] = verification(u1, u2, N_vec, do_plot)
%% VERIFICATION Spatial operator convergence via Method of Manufactured Solutions (MMS)
% Steps:
%   1. Compute exact convective and diffusive transport operators symbolically.
%   2. Loop through grid resolutions, evaluating analytical and numerical operators.
%   3. Locate maximum error positions and extract local staggered cell spacings (dx, dy).
%   4. Compute observed orders of accuracy p = log(e_prev/e_curr) / log(h_prev/h_curr).
%   5. Print formatted results table and display log-log convergence plot.



% Support both 1D row vector and 2D matrix (row 1: N, row 2: M)
if size(N_vec, 1) == 1
    N_vec = [N_vec; N_vec];
end
num_grids = size(N_vec, 2);

% 1. Exact transport operators (symbolic derivatives)
x = sym('x'); y = sym('y');
conv_sym_u = diff(u1^2, x) + diff(u1 * u2, y); % convective term for u
conv_sym_v = diff(u1 * u2, x) + diff(u2^2, y); % convective term for v
diff_sym_u = diff(u1, x, 2) + diff(u1, y, 2);  % diffusive term for u (Laplacian)
diff_sym_v = diff(u2, x, 2) + diff(u2, y, 2);  % diffusive term for v (Laplacian)

% Preallocate error and local step arrays
err_max_cu = zeros(1, num_grids);
err_max_cv = zeros(1, num_grids);
err_max_du = zeros(1, num_grids);
err_max_dv = zeros(1, num_grids);
h_cu = zeros(1, num_grids);
h_cv = zeros(1, num_grids);
h_du = zeros(1, num_grids);
h_dv = zeros(1, num_grids);

% 2. Grid refinement loop
for k = 1:num_grids
    N = N_vec(1, k);
    M = N_vec(2, k);
    mesh = create_mesh(N, M, 1, 1);
    
    % Node coordinates for u and v staggered grids
    nodes_u = cat(3, mesh.XU, mesh.YU);
    nodes_v = cat(3, mesh.XV, mesh.YV);

    % Analytical evaluations on staggered nodes
    ua = evaluation(u1, nodes_u);
    va = evaluation(u2, nodes_v);
    cu_analytic = evaluation(conv_sym_u, nodes_u);
    cv_analytic = evaluation(conv_sym_v, nodes_v);
    du_analytic = evaluation(diff_sym_u, nodes_u);
    dv_analytic = evaluation(diff_sym_v, nodes_v);
    
    % Discrete numerical operators
    cu = convective_u(ua, va, mesh);
    cv = convective_v(ua, va, mesh);
    du = diffusive_u(ua, mesh);
    dv = diffusive_v(va, mesh);

    % Interior error fields (excluding halo)
    e_cu = abs(cu_analytic(2:end-1, 2:end-1) - cu(2:end-1, 2:end-1));
    e_cv = abs(cv_analytic(2:end-1, 2:end-1) - cv(2:end-1, 2:end-1));
    e_du = abs(du_analytic(2:end-1, 2:end-1) - du(2:end-1, 2:end-1));
    e_dv = abs(dv_analytic(2:end-1, 2:end-1) - dv(2:end-1, 2:end-1));

    % Locate maximum errors and their node indices
    [err_max_cu(k), idx_cu] = max(e_cu(:));
    [err_max_cv(k), idx_cv] = max(e_cv(:));
    [err_max_du(k), idx_du] = max(e_du(:));
    [err_max_dv(k), idx_dv] = max(e_dv(:));

    [i_cu, ~] = ind2sub(size(e_cu), idx_cu); i_cu = i_cu + 1; % shift for halo
    [~, j_cv] = ind2sub(size(e_cv), idx_cv); j_cv = j_cv + 1;
    [i_du, ~] = ind2sub(size(e_du), idx_du); i_du = i_du + 1;
    [~, j_dv] = ind2sub(size(e_dv), idx_dv); j_dv = j_dv + 1;

    % Local staggered cell spacing at max error location (dx for u, dy for v)
    h_cu(k) = mesh.xp(i_cu + 1) - mesh.xp(i_cu);
    h_du(k) = mesh.xp(i_du + 1) - mesh.xp(i_du);
    h_cv(k) = mesh.yp(j_cv + 1) - mesh.yp(j_cv);
    h_dv(k) = mesh.yp(j_dv + 1) - mesh.yp(j_dv);
end

% 3. Observed convergence order
order_cu = nan(1, num_grids);
order_cv = nan(1, num_grids);
order_du = nan(1, num_grids);
order_dv = nan(1, num_grids);

for k = 2:num_grids
    order_cu(k) = log(err_max_cu(k-1) / err_max_cu(k)) / log(h_cu(k-1) / h_cu(k));
    order_cv(k) = log(err_max_cv(k-1) / err_max_cv(k)) / log(h_cv(k-1) / h_cv(k));
    order_du(k) = log(err_max_du(k-1) / err_max_du(k)) / log(h_du(k-1) / h_du(k));
    order_dv(k) = log(err_max_dv(k-1) / err_max_dv(k)) / log(h_dv(k-1) / h_dv(k));
end

% 4. Display formatted table
fprintf('\n===============================================================================================================\n');
fprintf('                                 MMS SPATIAL OPERATOR GRID CONVERGENCE RESULTS\n');
fprintf('===============================================================================================================\n');
fprintf('    Grid         hx         hy        ||e(Cu)||     Order    ||e(Cv)||     Order    ||e(Du)||     Order    ||e(Dv)||     Order\n');
fprintf('---------------------------------------------------------------------------------------------------------------\n');
for k = 1:num_grids
    grid_str = sprintf('%dx%d', N_vec(1, k), N_vec(2, k));
    if k == 1
        fprintf('%9s   %8.5f   %8.5f    %10.4e     ---     %10.4e     ---     %10.4e     ---     %10.4e     ---\n', ...
            grid_str, h_cu(k), h_cv(k), err_max_cu(k), err_max_cv(k), err_max_du(k), err_max_dv(k));
    else
        fprintf('%9s   %8.5f   %8.5f    %10.4e    %5.2f    %10.4e    %5.2f    %10.4e    %5.2f    %10.4e    %5.2f\n', ...
            grid_str, h_cu(k), h_cv(k), err_max_cu(k), order_cu(k), err_max_cv(k), order_cv(k), ...
            err_max_du(k), order_du(k), err_max_dv(k), order_dv(k));
    end
end
fprintf('===============================================================================================================\n\n');

% 5. Log-log convergence plot

h_fig = figure('Name', 'Spatial Grid Convergence Analysis of Operators', 'Color', 'w');

loglog(h_cu, err_max_cu, '-o', 'LineWidth', 2, 'DisplayName', 'Convective U'); hold on;
loglog(h_cv, err_max_cv, '-s', 'LineWidth', 2, 'DisplayName', 'Convective V');
loglog(h_du, err_max_du, '-^', 'LineWidth', 2, 'DisplayName', 'Diffusive U');
loglog(h_dv, err_max_dv, '-d', 'LineWidth', 2, 'DisplayName', 'Diffusive V');

% Reference slope: 2nd order
ref_h = h_du;
ref_2nd = err_max_du(1) * (ref_h ./ ref_h(1)).^2;
loglog(ref_h, ref_2nd, '-.k', 'LineWidth', 1.2, 'DisplayName', '\mathcal{O}(h^2)');

hold off;
grid on; box on;
xlabel('Local grid spacing h', 'FontSize', 12);
ylabel('Maximum Error ||e||_{\infty}', 'FontSize', 12);
title('Spatial Grid Convergence Analysis of Operators', 'FontSize', 14);
legend('Location', 'southeast', 'FontSize', 10);


% 6. Output struct
results.N_vec = N_vec;
results.h_cu = h_cu;
results.h_cv = h_cv;
results.h_du = h_du;
results.h_dv = h_dv;
results.err_max_cu = err_max_cu;
results.err_max_cv = err_max_cv;
results.err_max_du = err_max_du;
results.err_max_dv = err_max_dv;
results.order_cu = order_cu;
results.order_cv = order_cv;
results.order_du = order_du;
results.order_dv = order_dv;

end