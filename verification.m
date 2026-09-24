function [L2_error, e]=verification(mesh,f_analytic,f_numeric)

% Verify operators
% Verify diffusive operator

sy




[e,L2_error] = compute_error(mesh,f_analytic,f_numeric);

end

function [e,L2_error] = compute_error(mesh,f_analytic,f_numeric)
% We take the analytic function and evaluate it at all points where we have
% nodes of the mesh.
x_coord=mesh.x;
y_coord=mesh.y;
f_exact = f_analytic(x_coord, y_coord);
e = f_exact - f_numeric;
L2_error = sqrt(sum(e(:).^2) / numel(e));

end