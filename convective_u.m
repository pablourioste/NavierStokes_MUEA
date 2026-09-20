function cu = convective_u(u,v,hx,hy,mesh)
% cu = convective_u(u,v,hx,hy,mesh)
% Volume integral of the convective term of the x-momentum equation at
% every inner u node. Halo of cu is not updated.
% Author: Pablo Urioste // Marc Antich // Martí Esquerda

N = size(u, 1) - 2; % Número de celdas interiores en x
M = size(u, 2) - 2; % Número de celdas interiores en y

cu = zeros(size(u));

if nargin == 3
    L = hx;
    hx = L / N;
    hy = L / M;
end

use_mesh = (nargin >= 5 && isfield(mesh, 'xu') && isfield(mesh, 'yu'));

for i = 2:N+1
    if use_mesh
        dx_minus = mesh.xu(i)   - mesh.xu(i-1);
        dx_plus  = mesh.xu(i+1) - mesh.xu(i);
    else
        dx_minus = hx;
        dx_plus  = hx;
    end
    for j = 2:M+1
        if use_mesh
            dy = mesh.yu(j) - mesh.yu(j-1);
        else
            dy = hy;
        end

        ue = (u(i+1,j) + u(i,j)) / 2;
        uw = (u(i,j)   + u(i-1,j)) / 2;
        un = (u(i,j+1) + u(i,j)) / 2;
        us = (u(i,j)   + u(i,j-1)) / 2;

        Fe = (u(i+1,j)*dy + u(i,j)*dy) / 2;
        Fw = (u(i-1,j)*dy + u(i,j)*dy) / 2;
        Fn = (v(i,j)*dx_minus   + v(i+1,j)*dx_plus) / 2;
        Fs = (v(i,j-1)*dx_minus + v(i+1,j-1)*dx_plus) / 2;

        cu(i,j) = ue*Fe - uw*Fw + un*Fn - us*Fs;
    end
end
end
