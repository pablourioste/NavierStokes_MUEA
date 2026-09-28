function [delta_t, delta_tc, delta_td] = time_stability(mesh, u, v, nu, f)
%% TIME_STABILITY Evaluates CFL and diffusive stability limits for explicit time stepping.
% Reference: Topic 2, Slide 66 & Chapter 6 (Part C).
%
% Inputs:
%   mesh : staggered grid structure from create_mesh (contains N, M, xu, yv, etc.)
%   u    : u-velocity field (size (N+2)x(M+2), including halos)
%   v    : v-velocity field (size (N+2)x(M+2), including halos)
%   nu   : kinematic viscosity [m^2/s]
%   f    : (optional) safety factor, default = 0.1 (Slide 66)
%
% Outputs:
%   delta_t  : stable time step: delta_t = f * min(delta_tc, delta_td)
%   delta_tc : convective time limit (CFL condition)
%   delta_td : diffusive time limit (von Neumann condition)


N = mesh.N;
M = mesh.M;

% 1. Local cell dimensions for interior control volumes
% dx: column vector (N x 1)
% dy: column vector (M x 1)
dx = mesh.xu(2:N+1) - mesh.xu(1:N);
dy = mesh.yv(2:M+1) - mesh.yv(1:M);

% 2. Extract interior velocity magnitudes
u_int = abs(u(2:N+1, 2:M+1));
v_int = abs(v(2:N+1, 2:M+1));

% 3. Convective time limit (CFL)
% delta_tc <= min(dx / |u|) and min(dy / |v|) for each CV
term_cu = dx ./ u_int;
term_cv = (dy') ./ v_int;

delta_tc = min([term_cu(:); term_cv(:)]);

% 4. Diffusive time limit:
% General 2D formula accounting for rectangular/non-uniform aspect ratios:
term_du = dx.^2/nu;
term_dv = dy.^2/nu;

delta_td = 1/2* min([term_du(:); term_dv(:)]);
% 5. Overall stable time step with safety factor (Slide 66)
delta_t = f * min(delta_tc, delta_td);

end