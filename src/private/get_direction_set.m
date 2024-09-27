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
%   Another thing to note is that if the direction set is already orthogonal, we may change the order of the
%   directions in the direction set. This is because the permutation vector p in the QR factorization may change
%   the order of the directions. However, the directions in the direction set are still the same.
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
    else
        % Find a maximal linearly independent subset of the direction_set, and supplement this subset with
        % new directions to make a basis of the full space. The final direction_set will be this basis.
        % First, we use QR factorization with permutation to find a maximal linearly independent subset of the direction_set.
        % Then, we supplement this subset with some directions in Q to make a basis of the full space.
        % For details about how to implement QR factorization with permutation, see the following explanation.
        % Actually, QR factorization with permutation is implemented using Gram-Schmidt orthogonalization.
        % First, we select the longest vector in the direction_set as the first vector in the new direction set.
        % Then, we normalize this vector to be the first column of the Q matrix, denoted as q_1.
        % The second column of Q, denoted as q_2, is the unit vector along the longest vector in the remaining 
        % vectors {d_2, ..., d_m} after removing the projection of q_1 from them. Repeat this process to get
        % the Q matrix and the R matrix. According to the above process, we have the following properties in the R matrix:
        % 1. The diagonal elements of R are monotonically decreasing.
        % 2. {R_ii}^2 >= sum_{j = i}^{n} {R_kj}^2, where k = i+1, ..., m.
        [Q, R, p] = qr(direction_set, "vector");
        num_directions = size(direction_set, 2);
        is_independent = (abs(diag(R)) >= 1e-10);
        direction_set = [direction_set(:, p(is_independent)) Q(:, p(~is_independent)) Q(:, num_directions+1 : n)]; 
    end
end

% Finally, set D to [d_1, -d_1, ..., d_n, -d_n], where d_i is the i-th vector in direction_set.
D = NaN(n, 2*n);
D(:, 1:2:2*n-1) = direction_set;
D(:, 2:2:2*n) = -direction_set;

end

