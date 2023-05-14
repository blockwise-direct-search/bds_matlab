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

% Precondition: If debug_flag is true, then pre-conditions is operated on
% input. If input_correctness is false, then assert may let the code crash.
debug_flag = is_debugging();
if debug_flag
    precondition_divide_searching_set(m, nb);
end

num_direction_block = floor(m/nb); % Number of directions each block (average)
if mod(m, nb) == 0 % Every block will have the same length of indices.
    num_directions = num_direction_block*ones(nb,1);
else % 
    num_directions = [(num_direction_block+1)*ones(mod(m, nb), 1); num_direction_block*ones(nb-mod(m, nb), 1)];
end

% Use cell to avoid the situation that each block has different length 
index_searching_set = cell(1,nb);
for i = 1:nb
    index_searching_set(:,i) = {sum(num_directions(1:i-1))+1:1:sum(num_directions(1:i))};
end

% Postcondition: If debug_flag is true, then post-conditions is operated on
% output. If output_correctness is false, then assert will let code crash.
if debug_flag
    postcondition_divide_searching_set(index_searching_set, nb);
end

end
