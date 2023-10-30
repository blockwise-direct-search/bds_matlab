function D = rotate_directions(D, x_init, xval)
keyboard
if (norm(x_init-xval) > 0)
    [n, m] = size(D);
    D_base = D(:, 1:2:m-1);
    num_base = size(D_base, 2);
    D_new = zeros(size(D_base));

    % Arrange the inner products in a column vector manner.
    algebraic_sum_directions = vertcat((xval-x_init)'*D_base)';

    for i = 1:num_base
        D_new(:, i) = D_base(:, i:num_base)*algebraic_sum_directions(i:num_base);
    end
    keyboard
    % Perform Gram-Schmidt orthogonalization on the column vectors of the matrix.
    D = zeros(n, num_base);

    % Normalize the first column vector.
    D(:, 1) = D_new(:, 1) / norm(D_new(:, 1));

    % Orthogonalize the remaining column vectors.
    for i = 2:num_base
        d = D_new(:, i);
        % Subtract the projection of the previous vectors.
        d = d - D(:, 1:i-1) * (D(:, 1:i-1)' * D_new(:, i));
        % Normalize to obtain orthogonal vectors.
        D(:, i) = d / norm(d);
    end
    keyboard
    D_init = zeros(n, m);
    D_init(:, 1:2:m-1) = D;
    D_init(:, 2:2:m) = -D;
    D = D_init;
end

end

