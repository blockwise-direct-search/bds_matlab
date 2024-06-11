function D = get_direction_set(n, options)
%direction_set generates the set of polling directions.
%
%   If the user does not input OPTIONS or the OPTIONS does not contain the field of direction_set,
%   D will be [e_1, -e_1, ..., e_n, -e_n], where e_i is the i-th coordinate vector.
%   Otherwise, options.direction_set should be a n-by-n nonsingular matrix, and D will be
%   set to [d_1, -d_1, ..., d_n, -d_n], where d_i is the i-th column in options.direction_set.
%   If the columns of options.direction_set are almost linearly dependent, then we will revise direction_set 
%   in the following way.
%   1. Remove the directions whose norms are too small.
%   2. Find directions that are almost parallel. Then preserve the first one and remove the others.
%   3. Find a maximal linearly independent subset of the directions, and supplement this subset with
%      new directions to make a basis of the full space. The final direction set will be this basis. 
%      This is done by QR factorization.
%

% Set options to an empty structure if it is not provided.
if nargin < 2
    options = struct();
end

if ~isfield(options, "direction_set")

    % Set the default direction set.
    direction_set = eye(n);

else

    direction_set = options.direction_set;

    % Determine whether the direction set contains NaN or Inf values and replace those
    % elements with 0.
    if any(isnan(direction_set) | isinf(direction_set), "all")
        warning("Some directions contain NaN or inf, which are replaced with 0.");
        direction_set(isnan(direction_set) | isinf(direction_set)) = 0;
    end

    % Check whether the direction set is linearly independent by checking whether the smallest
    % eigenvalue of the direction set is larger than 1.0e-10. There may exists some cases
    % where some of the eigenvalues of the direction set are complex numbers. In this case, we compare
    % the norm of the smallest eigenvalue with 1.0e-10. Using eigs to calculate the smallest eigenvalue
    % may encounter some convergence issues. The warning message may be displayed like this:
    % "0 of the 1 requested eigenvalues converged. Eigenvalues that did not converge are NaN."
    % We set the maximum number of iterations to 1000 and the tolerance to eps to avoid these issues.
    if norm(eigs(direction_set, 1, 'smallestabs', 'MaxIter', 1000, 'Tolerance', eps)) < 1.0e-10
        error('The direction set is not linearly independent.');
    end

    % Remove the directions whose norms are too small.
    % We include the factor sqrt(n) in the following threshold to reflect the dimension of the
    % problem.
    shortest_direction_norm = 10*sqrt(n)*eps;
    direction_norms = vecnorm(direction_set);
    short_directions = (direction_norms < shortest_direction_norm);
    if any(short_directions)
        warning("The direction set contains directions shorter than %g. They are removed.", ...
            shortest_direction_norm);
        direction_set = direction_set(:, ~short_directions);
        direction_norms = direction_norms(~short_directions);
    end
    
    % Find those directions that are almost parallel. Preserve the first one appearing in the
    % direction set and remove the others.
    parallel_directions = (abs(direction_set'*direction_set) > (1 - 1.0e-10) * (direction_norms' * direction_norms));
    % Parallel_directions is a symmetric binary matrix, where the (i,j)-th element is 1 if and 
    % only if the i-th and j-th vectors are almost parallel in the direction_set. 
    % Triu(parallel_directions, 1) returns a matrix whose strict upper triangular part of 
    % parallel_directions and the other part consists of only zero. Applying `find` to this 
    % matrix, we can get the row and column indices of the 1's in the strict upper triangular
    % part of parallel_directions, the column indices being contained in parallel_direction_indices below. 
    % An index i appears in parallel_direction_indices if and only if there exists a j < i 
    % such that the i-th and j-th vectors are almost parallel in the direction_set. 
    % Hence we will remove the directions indexed by parallel_direction_indices.
    [~, parallel_direction_indices] = find(triu(parallel_directions, 1));
    % Remove the duplicate indices appearing in the parallel_direction_indices.
    % parallel_direction_indices = unique(parallel_direction_indices);
    % Directions not indexed by parallel_direction_indices are preserved.
    preserved_indices = ~ismember(1:size(direction_set, 2), parallel_direction_indices);
    direction_set = direction_set(:, preserved_indices);
    % If the direction set is empty, we set it to be the identity matrix. Indeed, the removal of 
    % almost parallel directions cannot make an nonempty direction set become empty. Thus the
    % following code can be moved above the removal of almost parallel directions. However, we 
    % keep it here to make the code more robust. 
    if isempty(direction_set)
        direction_set = eye(n);
    end
    
    % If rank(direction_set) < n, we add new columns to direction_set to make the rank become n.
    % We use QR factorization with permutation to find such columns. 

    % First, we need to subtract a maximum linearly independent set from direction_set by QR factorization.
    [~, R, P] = qr(direction_set);
    [row, col] = find(triu(P) == 1);
    % permuted_sigma records the permutation of the columns of direction_set. It maps the indices of
    % the columns of direction_set, which is the first column of permuted_sigma, to the indices of the
    % columns of the permuted direction_set, which is the second column of permuted_sigma.
    permuted_sigma = [row, col];
    % The linear independence of the columns of direction_set is equivalent to the upper triangular matrix
    % R. Since the diagonal elements of R are in decreasing order, we can know the length of the maximal 
    % linear independent system of direction_set by finding the first element in the diagonal
    % of R that is smaller than 1e-10. Then, by the permutation matrix P, we can know the original indices
    % of the columns of direction_set that are linearly independent.
    if ~isempty(1 : find(abs(diag(R)) < 1e-10, 1)) && length(1 : find(diag(R) < 1e-10, 1)) == 1
        R_truncate_index = find(diag(R) < 1e-10, 1);
        [~, AP_indices] = ismember(1:(R_truncate_index - 1), permuted_sigma(:, 2));
        direction_set = direction_set(:, sort(AP_indices));
    % else
    %     direction_set = eye(n);
    end
    
    % Note: Actually, the above code may influence the order of the columns of direction_set. However,
    % the order of the columns of direction_set does not matter since each block will be visited once in one iteration.
    % Thus, we can ignore the order of the columns of direction_set.
    
    % The following columns of Q will be added to direction_set.
    % 1. The corresponding diagonal elements of R are tiny.
    % 2. Columns m+1 to n with m being the number of columns in direction_set provided that m < n.
    [Q, R, ~] = qr(direction_set);
    [~, m] = size(direction_set);
    % deficient_columns contains the indices of the tiny diagonal elements of 
    % R(1:min(m, n), 1:min(m, n)). 
    % We must transpose diag(R) since vecnorm will return a row vector. Otherwise, the following
    % comparison will return a matrix due to the implicit expansion.
    deficient_columns = ~(abs(diag(R(1:min(m, n), 1:min(m, n)))))' > ...
        10*eps*max(m,n)*vecnorm(R(1:min(m,n), 1:min(m,n)));
    direction_set = [direction_set, Q(:, deficient_columns), Q(:, m+1:end)];

    % If rank(direction_set) < n, we add new columns to direction_set to make the rank become n.
    % We use QR factorization with permutation to find such columns. 
    % The following columns of Q will be added to direction_set.
    % 1. The corresponding diagonal elements of R are tiny.
    % 2. Columns m+1 to n with m being the number of columns in direction_set provided that m < n.
    [Q, R, ~] = qr(direction_set);
    [~, m] = size(direction_set);
    % deficient_columns contains the indices of the tiny diagonal elements of 
    % R(1:min(m, n), 1:min(m, n)). 
    % We must transpose diag(R) since vecnorm will return a row vector. Otherwise, the following
    % comparison will return a matrix due to the implicit expansion.
    deficient_columns = ~(abs(diag(R(1:m, 1:m)))) > 1.0e-10*vecnorm(R(1:m, 1:m))';
    direction_set = [direction_set, Q(:, deficient_columns), Q(:, m+1:end)];
end

% Finally, set D to [d_1, -d_1, ..., d_m, -d_m], where d_i is the i-th vector in direction_set.
[~, m] = size(direction_set);
D = NaN(n, 2*m);
D(:, 1:2:2*m-1) = direction_set;
D(:, 2:2:2*m) = -direction_set;

end

