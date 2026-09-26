function U =halo_update(u)

% Tanto N como M ya contempla que son +2
N = size(u,1); % Mida de u files
M = size(u,2); % % Mida de u columnes

U = size (N , M);

%  Iteraciones de dentro

for i=2:N-1
     for j=2:M-1
         U(i,j)=u(i,j);
    end
 end



% Update de las horizontales

for i=2:N-1
    U(i,1) = u(i,M-1);
    U(i,M) = u(i,2);
end

% Update de las horizontales

for j=2:M-1
    U(1,j) = u(N-1, j);
    U(N,j) = u(2, j);

end

% Update corner halo cells

U(1,1) = U(N-1, M-1);
U(1,M) = U(N-1, 2);
U(N,1) = U(2, M-1);
U(N,M) = U(2,2);

end