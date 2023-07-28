function index_searching_set = divide_searching_set(m, nb)
%DIVIDE_SEARCHING_SET Get indices in each block.
%   INDEX_SEARCHING_SET = DIVIDE_SEARCHING_SET(M, NB)
%   returns a cell that store indices in each block, where M directions and
%   NB blocks.
%
%   EXamples
%   M = 11; NB = 3;
%   Since 11/3 = 3...2,
%   then indices is [4 4 3].
%   Thus, INDEX_SEARCHING_SET is a cell, where
%   INDEX_SEARCHING_SET{1} = [1, 2, 3, 4];
%   INDEX_SEARCHING_SET{2} = [5, 6, 7, 8];
%   INDEX_SEARCHING_SET{3} = [9, 10, 11].
%

%   Preconditions: If debug_flag is true, then preconditions is to verify
%   input. If input_correctness is false, then assert may let the code crash.
debug_flag = is_debugging();
if debug_flag
    % Assert m is a positive integer.
    assert(isintegerscalar(m) && m>0);
    % Assert nb is a positive integer.
    assert(isintegerscalar(nb) && nb>0);
    % Assert the number of directions is greater than the number of
    % blocks. (Preprocess it)
    assert(nb<=m);
end

% Number of directions each block (average)
num_direction_block = floor(m/nb); 
% Every block will have the same length of indices.
if mod(m, nb) == 0 
    num_directions = num_direction_block*ones(nb,1);
else %
    num_directions = [(num_direction_block+1)*ones(mod(m, nb), 1); num_direction_block*ones(nb-mod(m, nb), 1)];
end

% Use cell to avoid the situation that each block has different length
index_searching_set = cell(1,nb);
for i = 1:nb
    index_searching_set(:,i) = {sum(num_directions(1:i-1))+1:1:sum(num_directions(1:i))};
end

% Postconditions: If debug_flag is true, then postconditions is to verify
% output. If output_correctness is false, then assert will let code crash.
if debug_flag
    assert(length(index_searching_set) == nb);
end

end
