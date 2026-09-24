function f_eval = evaluation(f_analytic, nodes)
    % f_analytic: Expresión o vector simbólico (ej: [u1; u2])
    % nodes: Matriz N x M x d (capa 1 es X, capa 2 es Y)
    % f_eval: Matriz N x M x num_componentes
    
    vars = symvar(f_analytic);
    
    % Extraer capas de coordenadas en celdas para pasar a matlabFunction
    n_vars = length(vars);
    args = cell(1, n_vars);
    for k = 1:n_vars
        args{k} = nodes(:, :, k);
    end
    
    % Evaluar cada componente del vector simbólico por separado
    num_comp = numel(f_analytic);
    [N, M, ~] = size(nodes);
    f_eval = zeros(N, M, num_comp);
    
    for c = 1:num_comp
        f_num = matlabFunction(f_analytic(c), 'Vars', vars);
        f_eval(:, :, c) = f_num(args{:});
    end
end