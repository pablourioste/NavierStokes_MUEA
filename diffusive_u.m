function du = diffusive_u(u,mesh)
% du = diffusive_u(u,hx,hy)
% Volume integral of the diffusive term of the x-momentum equation at
% every inner u node, without the viscosity (multiply by nu outside).
% Halo of du is not updated.
% Author: Pablo Urioste // Marc Antich // Martí Esquerda

du = zeros(size(u));
N=mesh.N;
M=mesh.M;

for i=2:1:N+1
    dx=mesh.xv(i+1)-mesh.xv(i);
    for j=2:1:M+1
        dy=mesh.yv(j)-mesh.yv(j-1);
        
        % Definition of velocities:
        uE=u(i+1,j);
        uP=u(i,j);
        uW=u(i-1,j);
        uN=u(i,j+1);
        uS=u(i,j-1);

        %Definition of distances:
        xP=mesh.xu(i);
        yP=mesh.yu(j);
        xE=mesh.xu(i+1);
        xW=mesh.xu(i-1);
        yN=mesh.yu(j+1);
        yS=mesh.yu(j-1);

        du_est=(uE-uP)/(xE-xP);
        du_west=(uP-uW)/(xP-xW);
        du_nord=(uN-uP)/(yN-yP);
        du_south=(uP-uS)/(yP-yS);

        du(i,j) = dy*(du_est-du_west) + dx*(du_nord-du_south);
        % Correction to eliminate the integration
        du(i,j) = du(i,j) / (dx*dy);
    end
end




end
