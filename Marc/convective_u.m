function cu = convective_u(u, v, N, M, L, H)

% Integral del término convectivo de u sobre el volumen de control
% staggered de u(i,j).
%   u, v : campos de velocidad (N+2)x(M+2), con el halo ya relleno
%   cu   : integral de volumen del convectivo (solo nodos interiores)


hx = L/N;
hy = H/M;

cu = zeros(N+2, M+2);

for i = 2:N+1
    for j = 2:M+1

        % Variable transportada en las caras (interpolación)
        ue = (u(i+1,j)+u(i,j))/2;
        uw = (u(i-1,j)+u(i,j))/2;
        un = (u(i,j+1)+u(i,j))/2;
        u_s = (u(i,j-1)+u(i,j))/2;

        v_n = (v(i+1,j)+v(i,j))/2;
        v_s = (v(i+1,j-1)+v(i,j-1))/2;

        % Flujos másicos por las caras (velocidad normal * longitud de la cara)
        Fe = ue*hy;
        Fw = uw*hy;
        Fn = v_n*hx;
        Fs = v_s*hx;

        cu(i,j) = ue*Fe - uw*Fw + un*Fn - u_s*Fs;

    end
end

end
