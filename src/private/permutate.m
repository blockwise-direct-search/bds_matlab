function [block_indices] = permutate(block_indices, options)
% If strategy is randomized_array, then permutate block_indices.

% Whether variable is an array.
assert(isnumeric(block_indices) && numel(block_indices) >= 1);
% Whether each of the elements of variable is a positive integer.
assert(all(block_indices(:) > 0) && all(mod(block_indices(:), 1) == 0));


if options.polling_blocks == "Randomized_array"
    nb = length(block_indices);
    block_indices = randperm(nb);
else 
    return;
end

% Whether output is an array.
assert(isnumeric(block_indices) && numel(block_indices) > 1);
% Whether each of the elements of output is a positive integer.
assert(all(block_indices(:) > 0) && all(mod(block_indices(:), 1) == 0));

end

