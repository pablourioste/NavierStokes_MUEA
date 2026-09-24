function cv = convective_v(u,v,mesh)
% cv = convective_v(u,v,hx,hy,mesh)
% Volume integral of the convective term of the y-momentum equation at
% every inner v node. Halo of cv is not updated.
% Author: Pablo Urioste // Marc Antich // Martí Esquerda

N = size(v, 1) - 2; % Número de celdas interiores en x
M = size(v, 2) - 2; % Número de celdas interiores en y

cv = zeros(size(v));

for i = 2:N+1
    dx = mesh.xu(i) - mesh.xu(i-1);
    for j = 2:M+1
        dy_plus  = mesh.yu(j+1) - mesh.yu(j);
        dy_minus = mesh.yu(j)   - mesh.yu(j-1);
        dy  = mesh.yp(j+1) - mesh.yp(j);   
       
        vn = (v(i,j+1) + v(i,j)) / 2;
        vs = (v(i,j-1) + v(i,j)) / 2;
        ve = (v(i+1,j) + v(i,j)) / 2;
        vw = (v(i-1,j) + v(i,j)) / 2;
        
        Fn = (v(i,j+1)*dx + v(i,j)*dx) / 2;
        Fs = (v(i,j-1)*dx + v(i,j)*dx) / 2;
        Fe = (u(i,j)*dy_minus   + u(i,j+1)*dy_plus) / 2;
        Fw = (u(i-1,j)*dy_minus + u(i-1,j+1)*dy_plus) / 2;
        
        cv(i,j) = vn*Fn - vs*Fs + ve*Fe - vw*Fw;
        % Correction to eliminate the integration
        cv(i,j) = cv(i,j) / (dx*dy);
    end
end
end


