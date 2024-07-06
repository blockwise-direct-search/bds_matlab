function p = crescent_10()

p.objective = @crescent;
p.x0 = [1, 1];

end



function f = crescent(x)

n = length(x);
f_index = zeros(2,1); 
for i = 1:2
    if i == 1
        for j = 1:n-1
            f_index(i) = f_index(i) + x(i)^2 + (x(i+1) - 1)^2 + x(i+1) - 1;
        end
    end
    if i == 2
        f_index(i) = 2;
    end
    f = f + -x(i)^2 - (x(i+1) - 1)^2 - x(i+1) + 1 + x(i+1) + 1;
end

end