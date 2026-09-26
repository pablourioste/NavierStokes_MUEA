function [Xu, Yu, Xv, Yv] = get_coord(N, M, L, H)

hx = L/N;
hy = H/M;

i = (1:N+2);
j = (1:M+2);

% Pression coordinates
Xp = (i-1.5)*hx;
Yp = (j-1.5)*hy;

% Vel u coordinates
Xu = Xp + 0.5*hx; % or Xu = (i-1)*hx;
Yu = Yp;

% Vel v coordinates
Xv = Xp;
Yv = Yp + 0.5*hy; % or Yv = (j-1)*hy;
end