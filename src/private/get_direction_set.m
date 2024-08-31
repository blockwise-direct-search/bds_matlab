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

    % Remove the directions whose norms are too small.
    % We include the factor sqrt(n) in the following threshold to reflect the dimension of the
    % problem.
    shortest_direction_norm = 10*sqrt(n)*eps;
    direction_norms = vecnorm(direction_set);
    short_directions = (direction_norms < shortest_direction_norm);
    % If the direction set is empty, the direction set should be set to an identity matrix. 
    if any(short_directions) && ~isempty(direction_set)
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

    % Todo: if the rank of direction set is already n, how does the qr factorization work?
    % we use QR factorization with permutation to make the rank of direction_set become n.
    % How does QR factorization with permutation work? It works under Gram-Schmidt orthogonalization.
    % Denote the columns of direction_set as [d_1, d_2, ..., d_m] where m is the number of columns of direction_set.
    % By Gram-Schmidt orthogonalization, we have 
    % d_i = (q_1^T d_1) q_1 + (q_2^T d_2) q_2 + ... + (q_i^T d_i) q_i
    %     = R_{1i} q_1 + R_{2i} q_2 + ... + R_{ii} q_i.
    % where q_i is the i-th column of Q and R is an upper triangular matrix.
    % So how to to make the diagonal elements of R decrease monotonically? First, we select the longest vector
    % as the first vector in Q, which is q_1. Then we select the longest vector after removing the projection
    % of q_1 from the original vectors as the second vector in Q, which is q_2. We continue this process until
    % we get the last vector q_m, which is the last column of Q. In this way, we can make the diagonal elements
    % of R decrease monotonically. Thus, we have R_{ii} > \sum_{j = i + 1}^{k} R_{j k} for i = 1, 2, ..., n, where k = i + 1, i + 2, ..., n.
    % If R_{ii} is tiny, then R_{jk} should be tiny for j = i + 1, i + 2, ..., n and k = i, i + 1, ..., n.
    % Therefore, the rank of direction_set is the number of non-tiny diagonal elements of R and the maximum
    % linearly independent subset of direction_set is [d_1, d_2, ..., d_r] where r is the last index of the
    % diagonal elements of R that are not tiny. In addition, [d_1, d_2, ..., d_r] can be spanned into
    % [d_1, d_2, ..., d_r, q_{r+1}, q_{r+2}, ..., q_n], which is a basis of the full space.
    % We need to use permutated vector p to reorder [d_1, ..., d_r] under the original order of direction_set.
    % We need to point out direction_set is not empty since we have set it to be the identity matrix if it is empty.
    % Thus, there exists at least one nonzero diagonal element in R to QR factorize direction_set.
    [Q, R, p] = qr(direction_set, "vector");
    is_independent = false(n, 1);
    is_independent(1:size(direction_set, 2)) = (abs(diag(R)) >= 1e-10);
    direction_set = [direction_set(:, p(is_independent)) Q(:, p(~is_independent))];

end

% Finally, set D to [d_1, -d_1, ..., d_m, -d_m], where d_i is the i-th vector in direction_set.
[~, m] = size(direction_set);
D = NaN(n, 2*m);
D(:, 1:2:2*m-1) = direction_set;
D(:, 2:2:2*m) = -direction_set;

end

