function dv = diffusive_v(v,mesh)
% dv = diffusive_v(v,mesh)
% Evaluates the diffusive term of the y-momentum equation at
% every inner v node (divided by cell volume), without viscosity (multiply by nu outside).
% Halo of dv is not updated.
% Author: Pablo Urioste // Marc Antich // Martí Esquerda // Iván Aguilar

dv = zeros(size(v));
N=mesh.N;
M=mesh.M;
for i=2:1:N+1
    dx=mesh.xu(i)-mesh.xu(i-1);
    for j=2:1:M+1
        dy=mesh.yu(j+1)-mesh.yu(j);

        % Velocity definitions:
        vE=v(i+1,j);
        vP=v(i,j);
        vW=v(i-1,j);
        vN=v(i,j+1);
        vS=v(i,j-1);

        % Coordinate distances:
        xP=mesh.xv(i);
        yP=mesh.yv(j);
        xE=mesh.xv(i+1);
        xW=mesh.xv(i-1);
        yN=mesh.yv(j+1);
        yS=mesh.yv(j-1);

        dv_est=(vE-vP)/(xE-xP);
        dv_west=(vP-vW)/(xP-xW);
        dv_nord=(vN-vP)/(yN-yP);
        dv_south=(vP-vS)/(yP-yS);

        dv(i,j) = dy*(dv_est-dv_west) + dx*(dv_nord-dv_south);
        % Divide by cell volume to obtain the differential operator
        dv(i,j) = dv(i,j) / (dx*dy);
    end
end

