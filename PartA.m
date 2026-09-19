%% Aerodinamica, Mecanica de Vol i Orbital - Project: Part A
% Author: Pablo Urioste // Marc Antich // Martí Esquerda
% Last version: 16/09/2026 --18:07

clear; clc; close all;

%% 1. Parámetros del Dominio
N = 3; % Número de celdas en x
M = 3; % Número de celdas en y
L = 1; % Longitud en x
H = 1; % Altura en y

%% 2. Generación y Visualización de Malla
mesh = create_mesh(N, M, L, H);

% Mostrar la malla en figura
plot_staggered_mesh(mesh);

%% 3. Definición de Funciones Analíticas de Velocidad
fu = @(x,y) x;
fv = @(x,y) x;

%% 4. Creación e Impresión del Campo de Prueba
F = build_test_field(N+2, M+2);
print_field(F, 'Campo de Prueba Original F');

%% 5. Actualización de Celdas Halo
% U = halo_update(F);
% print_field(U, 'Campo F con Halos Actualizados U');

%% 6. Evaluación del Campo de Velocidades
[u, v] = set_velocity_field(mesh.xu, mesh.yu, mesh.xv, mesh.yv, fu, fv);
u=halo_update(u);
v=halo_update(v);
print_field(u, 'Campo de Velocidad u');
print_field(v, 'Campo de Velocidad v');


function F = build_test_field(N,M)
% Creates an NxM matrix where F(i,j) encodes its own (i,j) position
F = zeros(N,M);
for i = 1:N
    for j = 1:M
        F(i,j) = i + j/10;
    end
end
end
