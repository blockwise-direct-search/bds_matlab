function D = get_searching_set(n, options)
   %SEARCHING_SET generates the set of polling directions.
   %
   %D = SEARCHING_SET(N) generates the set of polling directions
   %{e_1, -e_1, ..., e_n, -e_n}, which is canonical and represented in matrix form.
   %
   %D = SEARCHING_SET(N, OPTIONS) allows to provide options to the set
   %generation. Set OPTIONS.direction to "identity" to obtain
   %{e_1, ..., e_n, -e_1,..., -e_n}, represented in matrix form.
   %Set OPTIONS.direction to "canonical" to obtain {e_1, -e_1, ..., e_n, -e_n}, 
   %represented in matrix form.
   % 

% Set options to an empty structure if it is not provided.
if nargin < 2
    options = struct();
end

if isfield(options, "searching_set")

    % Remove the directions which their norm are too small from the searching set.
    searching_set = options.searching_set;
    shortest_direction_norm = 10*eps;
    direction_norms = sqrt(sum(searching_set.^2, 1));
    short_directions = (direction_norms < shortest_direction_norm);
    if any(short_directions)
        warning("The searching set contains directions shorter than %g. They are removed.", ...
            shortest_direction_norm);  
        searching_set = searching_set(:, ~short_directions);
    end
    
    % By QR factorization, the searching set will linearly span the full space.
    [Q, R] = qr(searching_set);
    [~, m] = size(Q);
    Q_indices = 1:m;
    Q_indices(abs(diag(R)) > 10*eps) = [];
    %searching_set = [searching_set Q(:, diag(R) <= eps)];
    searching_set = [searching_set, Q(:, Q_indices)]; 
    
    % Find those directions which are opposite or parallel in pairs. Preserve the
    % first one appearing in the searching set and remove the others.
    direction_norms = sqrt(sum(searching_set.^2, 1));
    G = searching_set'*searching_set; 
    opposite_directions = (G < -(1 - 1.0e-10) * (direction_norms' * direction_norms));   
    [~, repetition_directions_indices] = find(triu(opposite_directions));
    
    % Create a logical index vector to select the columns to preserve.
    preserved_indices = ~ismember(1:size(searching_set, 2), repetition_directions_indices);
    % Generate a new matrix by selecting the columns to preserve.
    D_clean = searching_set(:, preserved_indices);
    
    % Make the searching set positively span the full space.
    m = size(D_clean, 2);
    D = NaN(n, 2*m);
    D(:, 1:2:2*m-1) = D_clean;
    D(:, 2:2:2*m) = -D_clean;
else
    % Set the default searching set, which is identity.
    D = [eye(n) -eye(n)];
    % If the direction is not identity, then the searching set is canonical.
    if ~isfield(options, "direction") || options.direction ~= "identity"
        perm = zeros(2*n, 1);
        perm(1:2:2*n-1) = 1:n;
        perm(2:2:2*n) = n+1:2*n;
        D = D(:, perm);     
    end
end

end
