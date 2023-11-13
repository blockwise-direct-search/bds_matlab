function index_direction_set = divide_direction_set(m, nb)
%DIVIDE_DIRECTION_SET gets indices of the direction set in each block.
%   INDEX_direction_set = DIVIDE_DIRECTION_SET(M, NB) returns a cell that records indices of each 
%   block, where M is the number of directions and NB is the number of blocks. We try to make the
%   number of directions in each block as equal as possible. The last block may have fewer elements.
%
%   Example
%     M = 11, NB = 3.
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

% Calculate the member of directions of each block.
num_direction_block = floor(m/nb); 

if mod(m, nb) == 0 
    num_directions = num_direction_block*ones(nb,1);
else 
    % The last block may have fewer directions than others.
    num_directions = [(num_direction_block+1)*ones(mod(m, nb), 1); num_direction_block*ones(nb-mod(m, nb), 1)];
end

% Use cell instead of matrix in MATLAB to avoid the situation that each element has a different length.
index_direction_set = cell(1,nb);
for i = 1:nb
    index_direction_set(:,i) = {sum(num_directions(1:i-1))+1:1:sum(num_directions(1:i))};
end

% Detect whether the output is in the right type when debug_flag is true.
if debug_flag
    if length(index_direction_set) ~= nb
        error('The number of blocks of index_direction_set is not correct.');
    end
end

end
