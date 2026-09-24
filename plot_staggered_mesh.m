function plot_staggered_mesh(mesh)
% PLOT_STAGGERED_MESH Visualizes the staggered grid, CV walls, and nodes.
% Author: Pablo Urioste // Marc Antich // Martí Esquerda // Iván Aguilar

figure('Name', 'Staggered Grid Mesh Visualization', 'Color', 'w');
hold on;

hx = mesh.hx;
hy = mesh.hy;

% X and Y limits including halo layer
x_min = mesh.xp(1) - hx/2;
x_max = mesh.xp(end) + hx/2;
y_min = mesh.yp(1) - hy/2;
y_max = mesh.yp(end) + hy/2;

% 1. Plot Control Volume Boundaries (Cell Faces)
% Vertical cell faces (x = constant)
x_faces = mesh.xp - hx/2;
x_faces = [x_faces; mesh.xp(end) + hx/2];
for i = 1:length(x_faces)
    plot([x_faces(i), x_faces(i)], [y_min, y_max], '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.8);
end

% Horizontal cell faces (y = constant)
y_faces = mesh.yp - hy/2;
y_faces = [y_faces; mesh.yp(end) + hy/2];
for j = 1:length(y_faces)
    plot([x_min, x_max], [y_faces(j), y_faces(j)], '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.8);
end

% 2. Outline Physical Domain Boundary [0, L] x [0, H]
rectangle('Position', [0, 0, mesh.L, mesh.H], 'EdgeColor', 'k', ...
    'LineWidth', 2, 'LineStyle', '-');

% 3. Plot Nodes
% Pressure nodes (Cell Centers)
h_p = plot(mesh.XP(:), mesh.YP(:), 'o', 'MarkerEdgeColor', 'b', ...
    'MarkerFaceColor', [0.3 0.6 1], 'MarkerSize', 7, 'DisplayName', 'Pressure Nodes (P)');

% u-velocity nodes (x-faces)
h_u = plot(mesh.XU(:), mesh.YU(:), '>', 'MarkerEdgeColor', [0.8 0 0], ...
    'MarkerFaceColor', [1 0.3 0.3], 'MarkerSize', 6, 'DisplayName', 'u-velocity Nodes');

% v-velocity nodes (y-faces)
h_v = plot(mesh.XV(:), mesh.YV(:), '^', 'MarkerEdgeColor', [0 0.6 0], ...
    'MarkerFaceColor', [0.3 0.9 0.3], 'MarkerSize', 6, 'DisplayName', 'v-velocity Nodes');

% 4. Styling and Labels
axis equal;
xlim([x_min - 0.2*hx, x_max + 0.2*hx]);
ylim([y_min - 0.2*hy, y_max + 0.2*hy]);
xlabel('x');
ylabel('y');
title(sprintf('Arakawa-C Staggered Mesh (%dx%d Domain + Halo)', mesh.N, mesh.M));

legend([h_p, h_u, h_v], 'Location', 'northeastoutside');
box on;
hold off;
end