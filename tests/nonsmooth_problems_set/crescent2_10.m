function p = crescent2_10()

p.objective = @crescent2;
p.x0 = [1, 1];

end



function f = crescent2(x)

n = length(x);
f = 0;
for i = 1:n-1
    f = f + max((x(i)^2 + (x(i+1) - 1)^2 + x(i+1) - 1), (-x(i)^2 - (x(i+1) - 1)^2 - x(i+1) + 1 + x(i+1) + 1));
end

end