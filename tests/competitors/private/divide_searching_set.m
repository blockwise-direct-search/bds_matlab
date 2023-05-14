function index_searching_set = divide_searching_set(m, nb, cycling, polling)
%DIVIDE_SEARCHING_SET Get indices in each block.
%   INDEX_SEARCHING_SET = DIVIDE_SEARCHING_SET(M, NB, CYCLING, POLLING) 
%   returns a cell that store indices in each block, where M directions and
%   NB blocks. When CYCLING = 5 and POLLING = "opportunistic", the indices 
%   will be produced by the idea of sGauss-Seidel (symmetric Gauss Seidel).
% 
%   EXamples
%   M = 11; NB = 3; Either CYCLING ~= 5 or POLLING ~= "opportunistic"
%   Since 11/3 = 3...2,
%   then indices is [4 4 3].
%   Thus, INDEX_SEARCHING_SET is a cell, where
%   INDEX_SEARCHING_SET{1} = [1, 2, 3, 4];
%   INDEX_SEARCHING_SET{2} = [5, 6, 7, 8];
%   INDEX_SEARCHING_SET{3} = [9, 10, 11].
%
%   M = 11; NB = 3; CYCLING == 5 and POLLING == "opportunistic"
%   Since 11/3 = 3...2,
%   then indices is [4 4 3].
%   Thus, INDEX_SEARCHING_SET is a cell, where
%   INDEX_SEARCHING_SET{1} = [1, 2, 3, 4, 3, 2, 1];
%   INDEX_SEARCHING_SET{2} = [5, 6, 7, 8, 7, 6, 5];
%   INDEX_SEARCHING_SET{3} = [9, 10, 11, 10, 9].


% Precondition: If debug_flag is true, then pre-conditions is operated on
% input. If input_correctness is false, then assert may let the code crash.
debug_flag = is_debugging();
if debug_flag
    precondition_divide_searching_set(m, nb, cycling, polling);
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
    if cycling == 5 && ~strcmpi(polling, "complete") && num_directions(i)>1
        index_searching_set(:,i) = {[sum(num_directions(1:i-1))+1:1:sum(num_directions(1:i))...
            sum(num_directions(1:i))-1:-1:sum(num_directions(1:i-1))+1]};
    else
        index_searching_set(:,i) = {sum(num_directions(1:i-1))+1:1:sum(num_directions(1:i))};
    end
end

% Postcondition: If debug_flag is true, then post-conditions is operated on
% output. If output_correctness is false, then assert will let code crash.
if debug_flag
    postcondition_divide_searching_set(index_searching_set, nb);
end

end
