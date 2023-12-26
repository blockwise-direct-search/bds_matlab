function index_direction_set = divide_direction_set(n, num_blocks)
%DIVIDE_DIRECTION_SET gets indices of the direction set in each block.
%   INDEX_DIRECTION_SET = DIVIDE_DIRECTION_SET(n, num_blocks) returns a cell that records 
%   indices of directions in each block.
%   We try to make the number of directions in each block as even as possible. The last 
%   block may have one element fewer.
%
%   Example
%     n = 11, num_blocks = 3.
%     Since 11/3 = 3...2,
%     then indices are [4 4 3].
%     Thus, INDEX_DIRECTION_SET is a cell, where
%     index_direction_set{1} = [1, 12, 2, 13, 3, 14, 4, 15],
%     index_direction_set{2} = [5, 16, 6, 17, 7, 18, 8, 19],
%     index_direction_set{3} = [9, 20, 10, 21, 11, 22].
%

% Detect whether the input is given in the correct type when debug_flag is true.
debug_flag = is_debugging();
if debug_flag
    % n should be a positive integer.
    if ~isintegerscalar(n) || n <= 0    
        error('n is not a positive integer.');
    end
    % num_blocks should be a positive integer.
    if ~isintegerscalar(num_blocks) || num_blocks<=0    
        error('num_blocks is not a positive integer.');
    end
    % The number of blocks should not be greater than the number of variables.
    if n < num_blocks
        error('The number of blocks should not be greater than the number of variables.');
    end
end

% Calculate the number of directions of each block.
num_direction_each_block = floor(n/num_blocks); 

% num_direction_block indicates the number of directions of each block after averaging.
if mod(n, num_blocks) == 0 
    num_direction_block = num_direction_each_block*ones(num_blocks, 1);
else 
    % The last block may have fewer directions than others.
    num_direction_block = [(num_direction_each_block+1)*ones(mod(n, num_blocks), 1); num_direction_each_block*ones(num_blocks-mod(n, num_blocks), 1)];
end

% Use cell instead of matrix in MATLAB to avoid the number of directions in each block being different.
index_direction_set = cell(1, num_blocks);
for i = 1:num_blocks
    block_index_direction_set = sum(num_direction_block(1:i-1))+1:1:sum(num_direction_block(1:i));
    block_index_direction_set = 2*block_index_direction_set - 1;
    index_direction_set(:, i) = {reshape([block_index_direction_set; block_index_direction_set + 1], 1, [])};
end

% Detect whether the output is in the right type when debug_flag is true.
if debug_flag
    if length(index_direction_set) ~= num_blocks
        error('The number of blocks of index_direction_set is not correct.');
    end
end

end
