function p = cb32_10()

p.objective = @cb32;
p.x0 = [1, 1];

end



function f = cb32(x)

    n = length(x);
    f_index = zeros(3, 1);
    for i = 1:3
        if i == 1
            for j = 1:n-1
                f_index(i) = f_index(i) + x(j)^4 + x(j+1)^2;
            end
        end
        if i == 2
            for j = 1:n-1
                f_index(i) = f_index(i) + (2 - x(j))^2 + (2 - x(j+1))^2;
            end
        end
        if i == 3
            for j = 1:n-1
                f_index(i) = f_index(i) + 2*exp(-x(j) + x(j+1));
            end
        end
    end
    f = max(f_index);

end