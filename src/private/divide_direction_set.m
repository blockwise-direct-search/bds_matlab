function index_direction_set = divide_direction_set(num_directions, num_blocks)
%DIVIDE_DIRECTION_SET gets indices of the direction set in each block.
%   INDEX_DIRECTION_SET = DIVIDE_DIRECTION_SET(num_directions, num_blocks) returns a cell that records 
%   indices of each block, where num_directions	is the number of directions and num_blocks is the number
%   of directions. We try to make the number of directions in each block as even as possible. The last 
%   block may have one element fewer.
%
%   Example
%     num_directions = 11, num_blocks = 3.
%     Since 11/3 = 3...2,
%     then indices are [4 4 3].
%     Thus, INDEX_DIRECTION_SET is a cell, where
%     index_direction_set{1} = [1, 2, 3, 4],
%     index_direction_set{2} = [5, 6, 7, 8],
%     index_direction_set{3} = [9, 10, 11].
%

% Detect whether the input is given in the correct type when debug_flag is true.
debug_flag = is_debugging();
if debug_flag
    % num_directions should be a positive integer.
    if ~isintegerscalar(num_directions) || num_directions <= 0    
        error('num_directions is not a positive integer.');
    end
    % num_blocks should be a positive integer.
    if ~isintegerscalar(num_blocks) || num_blocks<=0    
        error('num_blocks is not a positive integer.');
    end
    % The number of directions should be greater than or equal to the number of blocks.
    if num_directions < num_blocks
        error('The number of directions should be greater than or equal to the number of blocks.');
    end
end

% Calculate the number of directions of each block.
num_direction_block = floor(num_directions/num_blocks); 

if mod(num_directions, num_blocks) == 0 
    num_directions = num_direction_block*ones(num_blocks, 1);
else 
    % The last block may have fewer directions than others.
    num_directions = [(num_direction_block+1)*ones(mod(num_directions, num_blocks), 1); num_direction_block*ones(num_blocks-mod(num_directions, num_blocks), 1)];
end

% Use cell instead of matrix in MATLAB to avoid the situation that each element has a different length.
index_direction_set = cell(1, num_blocks);
for i = 1:num_blocks
    index_direction_set(:, i) = {sum(num_directions(1:i-1))+1:1:sum(num_directions(1:i))};
end

% Detect whether the output is in the right type when debug_flag is true.
if debug_flag
    if length(index_direction_set) ~= num_blocks
        error('The number of blocks of index_direction_set is not correct.');
    end
end

end
