%% MMS diffusive
syms x y
L = 1; H = 1;

% Clacul funció simbolica
us =  cos(2*pi*x)*sin(2*pi*y);
vs = -2*sin(2*pi*x)*cos(2*pi*y); % per trencar simetria
% Cálcul del terme difusiu de l'equacio (part dreta eq momentum)
dus = diff(us,x,2) + diff(us,y,2);
dvs = diff(vs,x,2) + diff(vs,y,2);

% Fem estudi per diferents malles

Refinament  = [8 16 32 64 128];
err_u = zeros(size(Refinament));
err_v = zeros(size(Refinament));
h     = zeros(size(Refinament));

for k = 1:length(Refinament)
    N = Refinament(k);  M = N;
    hx = L/N;  hy = H/M;
    h(k) = hx;

   % 3) Sol exacte
    [du_ex, dv_ex] = set_vel_field(dus, dvs, N, M, L, H);

    % 4) Sol numerica
    [u, v] = set_vel_field(us, vs, N, M, L, H);
    du_num = diffusive_u(u, v, N, M, L, H) / (hx*hy);  
    dv_num = diffusive_v(u, v, N, M, L, H) / (hx*hy);

    % 5) Error: en els nodes interiors
    err_u (k) = max(max(abs( du_num(2:N+1,2:M+1) - du_ex(2:N+1,2:M+1) )));
    err_v (k) = max(max(abs( dv_num(2:N+1,2:M+1) - dv_ex(2:N+1,2:M+1) )));
end


figure
loglog(h, err_u, 'o-', h, err_v, 's--')
hold on
loglog(h, err_u(1)*(h/h(1)).^2, 'k:')   % referencia pendiente 2
xlabel('h [m]'); ylabel('max |error|')
legend('diffusive u', 'diffusive v', 'pendiente 2')
grid on

p_u = polyfit(log(h), log(err_u), 1);  orden_u = p_u(1);
p_v = polyfit(log(h), log(err_v), 1);  orden_v = p_v(1);
