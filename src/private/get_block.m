function block_index = get_block(k, nb, hist, polling_outer)

hist_block = hist.block;

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
    case {"Randomized"}
            block_index = round(rand(1,nb-1)) + 1;
end

end

