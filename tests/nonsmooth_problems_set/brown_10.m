function p = brown_10()

p.objective = @brown;
p.x0 = [1, 1];

end



function f = brown(x)

n = length(x);
f = 0;
for i = 1:n-1
    f = f + abs(x(i))^(x(i+1)^2 + 1) + abs(x(i+1))^(x(i)^2 + 1);
end

end