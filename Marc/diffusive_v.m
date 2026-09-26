function dv = diffussive_v(u, v, N, M, L, H)

% Integral del término difusivo de v sobre el volumen de control
% staggered de u(i,j).
%   u, v : campos de velocidad (N+2)x(M+2), con el halo ya relleno
%   dv   : integral de volumen del difusivo (solo nodos interiores)


hx = L/N;
hy = H/M;

dv = zeros(N+2, M+2);

for i = 2:N+1
    for j = 2:M+1
        dvdux_e = (v(i+1,j)-v(i,j))/hx;
        dvdx_w = (v(i,j)-v(i-1,j))/hx;
        dvdy_n = (v(i,j+1)-v(i,j))/hy;
        dvdy_s = (v(i,j)-v(i,j-1))/hy;

        dv(i,j) = dvdux_e*hy-dvdx_w*hy+dvdy_n*hx-dvdy_s*hx;

    end
end

end