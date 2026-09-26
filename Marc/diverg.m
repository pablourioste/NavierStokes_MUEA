function d = diverg (u, v, N, M, L, H)

hx = L/N;
hy = H/M;

d = zeros(N+2, M+2);

for i=2:N+1
    for j=2:M+1
        d(i,j)= (u(i,j)-u(i-1,j))*hy + (v(i,j)-v(i,j-1))*hx;
    end

end

end
