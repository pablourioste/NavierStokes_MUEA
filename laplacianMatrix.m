function A = laplacianMatrix(mesh)
%% LAPLACIANMATRIX Constructs the discrete Laplacian matrix with periodic BCs.
% Author: Pablo Urioste // Marc Antich // Martí Esquerda // Iván Aguilar
% Supports non-uniform rectangular grids directly from the mesh structure.

N = mesh.N;
M = mesh.M;

A = sparse(N*M, N*M);


for j = 2:M+1
    for i = 2:N+1
        % 1D index of current cell (equation row)
        p = (j - 2)*N + (i - 1);
        
        % East neighbor (periodic wrap-around at right boundary)
        if i < N+1
            east = p + 1;
        else
            east = p - (N - 1);
        end
        
        % West neighbor (periodic wrap-around at left boundary)
        if i > 2
            west = p - 1;
        else
            west = p + (N - 1);
        end
        
        % North neighbor (periodic wrap-around at top boundary)
        if j < M+1
            north = p + N;
        else
            north = p - (M - 1)*N;
        end
        
        % South neighbor (periodic wrap-around at bottom boundary)
        if j > 2
            south = p - N;
        else
            south = p + (M - 1)*N;
        end
        
        % Flux weights for finite volume Laplacian
    
        % Cell face dimensions (surface areas)
        dx = mesh.xu(i) - mesh.xu(i-1);
        dy = mesh.yv(j) - mesh.yv(j-1);
        
        % Distances between neighboring cell centers
        dx_e = mesh.xp(i+1) - mesh.xp(i);
        dx_w = mesh.xp(i)   - mesh.xp(i-1);
        dy_n = mesh.yp(j+1) - mesh.yp(j);
        dy_s = mesh.yp(j)   - mesh.yp(j-1);
        
        cx_e = dy / dx_e;
        cx_w = dy / dx_w;
        cy_n = dx / dy_n;
        cy_s = dx / dy_s;
        cp   = -(cx_e + cx_w + cy_n + cy_s);
        
        % Populate row p of matrix A
        A(p, p)     = cp;
        A(p, east)  = cx_e;
        A(p, west)  = cx_w;
        A(p, north) = cy_n;
        A(p, south) = cy_s;
    end
end

% Remove singularity for periodic Poisson equation (slide 56)
A(1, 1) = A(1, 1) - 1;
end