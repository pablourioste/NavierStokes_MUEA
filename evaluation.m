function f_eval = evaluation(f_analytic, nodes)
    % EVALUATION Evaluates symbolic expression or vector over grid nodes.
    % Input:
    %   f_analytic - Symbolic expression or vector (e.g., [u1; u2])
    %   nodes      - Array of size N x M x d (layer 1: X, layer 2: Y)
    % Output:
    %   f_eval     - Evaluated numerical array of size N x M x num_components
    % Author: Pablo Urioste // Marc Antich // Martí Esquerda // Iván Aguilar
    
    vars = symvar(f_analytic);
    
    % Extract coordinate layers into cells for matlabFunction input
    n_vars = length(vars);
    args = cell(1, n_vars);
    for k = 1:n_vars
        args{k} = nodes(:, :, k);
    end
    
    % Evaluate each component of the symbolic vector individually
    num_comp = numel(f_analytic);
    [N, M, ~] = size(nodes);
    f_eval = zeros(N, M, num_comp);
    
    for c = 1:num_comp
        f_num = matlabFunction(f_analytic(c), 'Vars', vars);
        f_eval(:, :, c) = f_num(args{:});
    end
end