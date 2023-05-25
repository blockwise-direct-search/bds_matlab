function block_index = get_block(nb, hist, blocks_strategy, block_indices)

% Extreme case: How about nb = 1? 

block_array_init = 1:nb;

block_hist = hist.block;
if nb ~= 1
    
    % Initialize hist_block as a zero vector. Thus, number of
    % nonzero elements will be number of blocks that having been
    % visited
    nb_visited = nnz(block_hist);
    
    switch blocks_strategy
        % Gauss-Seidel
        case {"Gauss-Seidel"}
            block_index = mod(nb_visited, nb)+1;
        
        % symmetric Gauss-Seidel    
        case {"sGauss-Seidel"}
            if nb_visited == 0 || nb_visited == 1
                block_index = nb_visited+1;
            else
                block_index_last = find(block_hist, 1, 'last');
                block_last = block_hist(block_index_last);
                block_second_last = block_hist(block_index_last-1);
                if block_last == 1 || (block_second_last < block_last && block_last ~= nb)
                    block_index = block_last+1;
                else
                    block_index = block_last-1;
                end
            end
            
        % Randomly produce an index from 1:nb.
        case {"Randomized_block_index"}
            block_index = round(rand(1,1)*(nb-1)+1);
        
        % Randomly produce an index from 1:nb, which is different from the 
        % last index having been visited
        case {"Randomized_no_repetition"}            
            if nb_visited == 0
                block_index = 1;
            else
                block_index_last = find(block_hist, 1, 'last');
                block_last = block_hist(block_index_last);
                block_array_init(block_last) = [];
                index = round(rand(1,1)*(nb-2)+1);
                block_index = block_array_init(index);
            end
            
        % Randomly produce an array from 1:nb. The following nb blocks going 
        % to be visited is the array.
        case {"Randomized_block_array"}
            index = mod(nb_visited, nb);
            block_index = block_indices(index)+1;

    end
else
   block_index = 1; 
end

end

