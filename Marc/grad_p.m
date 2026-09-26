function [gx, gy] = grad_p(p, N, M, L, H)
% [gx, gy] = grad_p(p, N, M, L, H)
% Gradient of a centred scalar field, evaluated at the staggered nodes.
%   p    : centred scalar field (N+2)x(M+2), halo already updated
%   N, M : number of interior cells in x and y
%   L, H : domain size in x and y [m]
%   gx   : dp/dx at u-nodes (N+2)x(M+2), interior only
%   gy   : dp/dy at v-nodes (N+2)x(M+2), interior only
% Point values (NOT integrated over the CV): used directly in u = up - grad(p)

hx = L/N;
hy = H/M;

gx = zeros(N+2, M+2);
gy = zeros(N+2, M+2);

for i = 2:N+1
    for j = 2:M+1
        gx(i,j) = (p(i+1,j)-p(i,j))/hx;      % dp/dx at u(i,j)
        gy(i,j) = (p(i,j+1)-p(i,j))/hy ;      % dp/dy at v(i,j)
    end
end

end


% ((p(i+1,j)-p(i,j))*hy)/hx + ((p(i-1,j)-p(i,j))*hy)/hx +
% ((p(i,j+1)-p(i,j))*hx)/hy + ((p(i,j-1)-p(i,j))*hx)/hy;Para el LAPLACIANO