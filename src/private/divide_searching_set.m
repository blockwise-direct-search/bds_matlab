function index_searching_set = divide_searching_set(m, nb)
%DIVIDE_SEARCHING_SET Get indices of the searching set in each block.
%   INDEX_SEARCHING_SET = DIVIDE_SEARCHING_SET(M, NB)
%   returns a cell that record indices of each block, where M is the number of directions
%   and the number of blocks is NB.
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

%   Preconditions: If debug_flag is true, then preconditions are to verify
%   input. If input_correctness is false, then assert may let the code crash.
debug_flag = is_debugging();
if debug_flag
    % m should be a positive integer.
    if ~isintegerscalar(m) || m<=0    
        error('m is not a positive integer.');
    end
    % nb should be a positive integer.
    if ~isintegerscalar(nb) || nb<=0    
        error('nb is not a positive integer.');
    end
    % The number of directions should be greater than or equal to the number of blocks.
    if m<nb
        error('The number of directions should be greater than or equal to the number of blocks.');
    end
end

% Calculate the mumber of directions of each block (average roughly).
num_direction_block = floor(m/nb); 

if mod(m, nb) == 0 
    num_directions = num_direction_block*ones(nb,1);
else 
    % The last block may have less directions than others.
    num_directions = [(num_direction_block+1)*ones(mod(m, nb), 1); num_direction_block*ones(nb-mod(m, nb), 1)];
end

% Use cell instead of matrix in MATLAB to avoid the situation that each block has different length.
index_searching_set = cell(1,nb);
for i = 1:nb
    index_searching_set(:,i) = {sum(num_directions(1:i-1))+1:1:sum(num_directions(1:i))};
end

% Postconditions: If debug_flag is true, then postconditions are to verify
% output. If output_correctness is false, then assert will let code crash.
if debug_flag
    if length(index_searching_set) ~= nb
        error('The number of blocks of index_searching_set is not correct.');
    end
end

end
