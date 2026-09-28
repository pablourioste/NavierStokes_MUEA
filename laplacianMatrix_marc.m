function A = laplacianMatrix_marc(mesh)
%% LAPLACIANMATRIX_MARC Constructs the discrete Laplacian matrix with periodic BCs.
% Author: Pablo Urioste // Marc Antich // Martí Esquerda // Iván Aguilar
% Computes neighbor indices using Marc's index map method (vector2field + halo_update).
% Supports non-uniform rectangular grids directly from the mesh structure.

N = mesh.N;
M = mesh.M;

% 1) Index map: K(i,j) = position of cell (i,j) in the algebraic vector.
%    vector2field places the numbers 1..N*M in the interior cells and
%    halo_update automatically wraps the periodic boundary neighbours.
K = vector2field((1:N*M)', mesh);
K = halo_update(K);

% 2) Sparse matrix initialization
A = sparse(N*M, N*M);

for j = 2:M+1
    for i = 2:N+1
        % Cell and neighbor indices from index map K (Marc's method)
        p     = K(i, j);
        east  = K(i+1, j);
        west  = K(i-1, j);
        north = K(i, j+1);
        south = K(i, j-1);
        
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
