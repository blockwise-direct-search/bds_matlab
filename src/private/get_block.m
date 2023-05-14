function block_index = get_block(k, nb, hist, polling_outer)

% Extreme case: How about nb = 1? 

hist_block = hist.block;
if nb ~= 1
    switch polling_outer
        case {"Gauss-Seidel"}
            % Gauss-Seidel
            if hist_block(length(hist_block)) == 0 || hist_block(k) == nb
                block_index = 1;
            else
                block_index = hist_block(k)+1;
            end
        case {"sGauss-Seidel"}
            % symmetric Gauss-Seidel
            if hist_block(length(hist_block)) == 0 || hist_block(k) == 1 || ...
                    (hist_block(k-1)< hist_block(k) && hist_block(k) ~= nb)
                block_index = hist_block(k)+1;
            else
                block_index = hist_block(k)-1;
            end
        case {"Randomized_block_index"}
            block_index = round(rand(1,1)*(nb-1)+1);
        case {"Randomized_no_repetition"}
            for i = 1:10000
                block_index = round(rand(1,1)*(nb-1)+1);
                if block_index ~= hist_block(k)
                    break;
                end
            end
    end
else
   block_index = 1; 
end

end

