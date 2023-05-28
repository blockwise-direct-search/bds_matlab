function [block_indices] = permutate(block_indices, options)
% If strategy is randomized_array, then permutate block_indices.


% Preconditions: If debug_flag is true, then pre-conditions is operated on
% input. If input_correctness is false, then assert may let the code crash.
debug_flag = is_debugging();

if debug_flag
    % Whether variable is an array.
    assert(isnumeric(block_indices) && numel(block_indices) >= 1);
    % Whether each of the elements of variable is a positive integer.
    assert(all(block_indices(:) > 0) && all(mod(block_indices(:), 1) == 0));
end

if options.blocks_strategy == "Randomized_array"
    nb = length(block_indices);
    block_indices = randperm(nb);
else 
    return;
end

% Postconditions: If debug_flag is true, then post-conditions is operated on
% output. If output_correctness is false, then assert will let code crash.
if debug_flag
    % Whether output is an array.
    assert(isnumeric(block_indices) && numel(block_indices) >= 1);
    % Whether each of the elements of output is a positive integer.
    assert(all(block_indices(:) > 0) && all(mod(block_indices(:), 1) == 0));
end

end

