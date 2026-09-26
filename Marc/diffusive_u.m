function du = diffussive_u(u, v, N, M, L, H)

% Integral del término difusivo de u sobre el volumen de control
% staggered de u(i,j).
%   u, v : campos de velocidad (N+2)x(M+2), con el halo ya relleno
%   du   : integral de volumen del difusivo (solo nodos interiores)


hx = L/N;
hy = H/M;

du = zeros(N+2, M+2);

for i = 2:N+1
    for j = 2:M+1
        dudux_e = (u(i+1,j)-u(i,j))/hx;
        dudx_w = (u(i,j)-u(i-1,j))/hx;
        dudy_n = (u(i,j+1)-u(i,j))/hy;
        dudy_s = (u(i,j)-u(i,j-1))/hy;

        du(i,j) = dudux_e*hy-dudx_w*hy+dudy_n*hx-dudy_s*hx;

    end
end

end
