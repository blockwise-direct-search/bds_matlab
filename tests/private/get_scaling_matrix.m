function sigma = get_scaling_matrix(n,options)
% Produce an scaling_matrix (diagonal)
% factor = options.scaling_matrix_factor;
sigma = eye(n);
factor = options.scaling_matrix_factor;
for i = 1:n
    sigma(i,i) = factor^i;
end

end

