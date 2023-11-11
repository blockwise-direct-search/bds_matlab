function D = get_direction_set(n, options)
%direction_set generates the set of polling directions.
%
%   If the user does not input OPTIONS or the OPTIONS does not contain direction_set.
%   D will be [e_1, -e_1, ..., e_n, -e_n], where e_i is the i-th coordinate vector.
%   Otherwise, options.direction_set should be a matrix of n rows, and set D to
%   [d_1, -d_1, ..., d_m, -d_m], where d_i is the i-th column of direction_set, 
%   and m is the number of columns of direction_set. Before doing so, we revise 
%   direction_set in the following way.
%   1. Remove the directions whose norms are too small. 
%   2. Find directions that are almost parallel. Then preserve the first one and remove the others.
%   3. Make the direction set linearly span the full space by adding new columns if necessary.
%      This is done by QR factorization.
%

% Set options to an empty structure if it is not provided.
if nargin < 2
    options = struct();
end

if ~isfield(options, "direction_set")

    % Set the default direction set.
    D = NaN(n, 2*n);
    D(:,1:2:2*n-1) = eye(n);
    D(:,2:2:2*n) = -eye(n);

else

    direction_set = options.direction_set;
    % Determine whether the direction set contains NaN or Inf values and replace those
    % elements with 0.
    if any(isnan(direction_set) | isinf(direction_set), "all")
        warning("Some directions contain NaN or inf. They are replaced with 0.");
        direction_set(isnan(direction_set) | isinf(direction_set)) = 0;
    end

    % Remove the directions whose norms are too small.
    % We include the factor sqrt(n) in the following threshold to reflect the dimension of the problem.
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
    % Parallel_directions is a symmetric matrix, whose diagonal is 1. Triu(parallel_directions, 1) returns
    % the indices of the upper triangular part of the matrix parallel_directions, setting the diagonal as zero.
    % By using find, we can get the row and column indices of the elements in the upper triangular part of
    % parallel_directions in pair, whose boolean value is 1. The column indices are the indices of the parallel
    % directions that do not appear for the first time in the direction set.
    [~, parallel_directions_indices] = find(triu(parallel_directions, 1));
    % Remove the duplicate indices appearing in the parallel_directions_indices.
    parallel_directions_indices = unique(parallel_directions_indices);
    
    % We remove the directions that are almost parallel, not appearing for the first time.
    preserved_indices = ~ismember(1:size(direction_set, 2), parallel_directions_indices);
    % Keep the remaining directions.
    direction_set = direction_set(:, preserved_indices);
    % If the direction set is empty, we set it to be the identity matrix.
    if isempty(direction_set)
        direction_set = eye(n);
    end

    % If rank(direction_set) < n, we add new columns to direction_set to make the rank become n.
    % We use QR factorization with permutation to find such columns. 
    [Q, R, ~] = qr(direction_set);
    [~, m] = size(R);
    
    % We add some columns from Q in direction set, where the corresponding index of diagonal elements of R 
    % are relatively small. When we determine which elements are too small, consider two types of matrix R 
    % separately.
    % 1. R is a vector. Then we only check the absolute value of the first element of R. 
    % 2. R is a matrix but not vector. Then we consider the absolute value of the diagonal elements of R.
    % we transpose diag(R) since the vecnorm will return a column vector.
    if min(m, n) == 1
        rank_direction_set_clean = sum(abs(R(1)) > 10*eps*max(m,n)*norm(R));
    else
        % rank_direction_set_clean = sum(abs(diag(R))' > 10*eps*max(m,n)*vecnorm(R(1:min(m,n), 1:min(m,n))));
        rank_direction_set_clean = find((abs(diag(R))' > 10*eps*max(m,n)*vecnorm(R(1:min(m,n), 1:min(m,n)))) == 0);
    end
    % To make D linearly span the full space, we add some columns in Q if necessary. 
    direction_set = [direction_set, Q(:, rank_direction_set_clean)];

    % Make the direction set positively span the full space.
    [direction_set_row_num, direction_set_col_num] = size(direction_set);
    D = NaN(direction_set_row_num, 2*direction_set_col_num);
    D(:,1:2:size(D, 2)-1) = direction_set;
    D(:,2:2:size(D, 2)) = -direction_set;

end

end
