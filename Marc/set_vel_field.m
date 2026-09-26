function [u,v] =  set_vel_field (us, vs, N, M, L, H)

% Em passo la funcio syms a normal (MatlabFunction)
syms x y
uf = matlabFunction(us, 'Vars', [x y]);
vf = matlabFunction(vs, 'Vars', [x y]);

[Xu, Yu, Xv, Yv] = get_coord(N, M, L, H);


u=zeros(N+2, M+2);
v=zeros(N+2, M+2);

for i=1:N+2
    for j=1:M+2
        u(i,j)=uf(Xu(i) , Yu(j));
        v(i,j)=vf(Xv(i) , Yv(j));
    end
end


% u(2:N+1, 2:M+1) = uf(Xu, Yu);
% v(2:N+1, 2:M+1) = vf(Xv, Yv);
% 
% u([1 end], :) = u([2 end-1], :);
% v(:, [1 end]) = v(:, [2 end-1]);



end