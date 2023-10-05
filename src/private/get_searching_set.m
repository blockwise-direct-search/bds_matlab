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
