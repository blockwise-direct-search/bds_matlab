function D = get_searching_set(n, options)
   %SEARCHING_SET generates the set of polling directions.
   %
   %D = SEARCHING_SET(N) generates the set of polling directions
   %{e_1, -e_1, ..., e_n, -e_n}, which is canonical and represented in matrix form.
   %
   %D = SEARCHING_SET(N, OPTIONS) allows to provide options to the set
   %generation. Set OPTIONS.direction to "identity" to obtain
   %{e_1, ..., e_n, -e_1,..., -e_n}, represented in matrix form.
   %

% Set options to an empty structure if it is not provided.
if nargin < 2
    options = struct();
end

if isfield(options, "searching_set")
    searching_set = options.searching_set;
    if rank(searching_set) < n
        warning("The rank of the searching set is advised to be n.");
    else
        if size(searching_set, 2) == n
            D = NaN(n, 2*n);
            D(:, 1:2:2*n-1) = searching_set;
            D(:, 2:2:2*n) = -searching_set;
        else
            [n, m] = size(searching_set); 
            % Generate indices for all possible combinations of vectors.
            [indices1, indices2] = meshgrid(1:m, 1:m);
            combinations = [indices1(:), indices2(:)];
            % Find the indices of vector combinations that sum to zero.
            zero_indices = all(searching_set(:, combinations(:, 1)) + ...
                searching_set(:, combinations(:, 2)) == 0, 1);
            % Extract the indices of vector combinations that sum to zero.
            index_matrix = combinations(zero_indices, :);
            % Remove duplicate index pairs.
            index_matrix = unique(sort(index_matrix, 2), 'rows');
            index_matrix_row = reshape(index_matrix.', 1, []);
            D_pair = searching_set(:, index_matrix_row);
            % Obtain the column indices which are not included into the 
            % vector combinations that sum to zero.
            all_indices = 1:m;
            preserve_indices = setdiff(all_indices, unique(index_matrix(:)));
            if ~isempty(preserve_indices)
                preserve_indices_length = length(preserve_indices);
                D_preserve = NaN(n, 2*preserve_indices_length);
                D_preserve_initial = searching_set(:, preserve_indices);
                D_preserve(:, 1:2:2*preserve_indices_length-1) = D_preserve_initial;
                D_preserve(:, 2:2:2*preserve_indices_length) = -D_preserve_initial;
                D = [D_pair D_preserve];
            else
                D = D_pair;
            end
        end
    end
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
