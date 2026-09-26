syms x y
L = 1; H = 1;
N = 32; M = 32;
hx = L/N; hy=H/M;

% 1) Camp
us =  cos(2*pi*x)*sin(2*pi*y);
vs = -sin(2*pi*x)*cos(2*pi*y);

% 2) Convectivo EXACTO (simbólico, forma conservativa)
cus = diff(us*us, x) + diff(us*vs, y);
cvs = diff(vs*us, x) + diff(vs*vs, y);

% 3) Exacto evaluado en los nodos (u-nodos para cus, v-nodos para cvs)
% Solución exacta de lo que s eva a comparar
[cu_ex, cv_ex] = set_vel_field(cus, cvs, N, M, L, H);

% 4) Numérico (aproximació)
[u, v] = set_vel_field(us, vs, N, M, L, H);
cu_num = convective_u(u, v, N, M, L, H) / (hx*hy);   
cv_num = convective_v(u, v, N, M, L, H) / (hx*hy);

% 5) Error: solo nodos interiores
err_u = max(max(abs( cu_num(2:N+1,2:M+1) - cu_ex(2:N+1,2:M+1) )))
err_V = max(max(abs( cv_num(2:N+1,2:M+1) - cv_ex(2:N+1,2:M+1) )))
