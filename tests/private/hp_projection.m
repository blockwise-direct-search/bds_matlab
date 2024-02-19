function x_projected = hp_projection(x)
    % TODO: IS THIS A PROJECTION?
    x_projected = x(:);
    x_projected(1) = max(x_projected(1), 1);
    x_projected(2) = max(min(x_projected(2), 1 - eps), eps);
    x_projected(3) = max(x_projected(3), 0);
    x_projected(4) = max(max(x_projected(4), eps), x_projected(3));
    x_projected(5) = max(x_projected(5), x_projected(4));
    return;
end
