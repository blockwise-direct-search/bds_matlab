function p = cb3_10()

p.objective = @cb3;
p.x0 = [1, 1];

end



function f = cb3(x)

n = length(x);
f_index = zeros(n-1, 1);
for i = 1:n-1
    f_index(i) = max(max(x(i)^4 + x(i+1)^4, (2 - x(i))^2 + (2 - x(i+1))^2), 2*exp(-x(i) + x(i+1)));
end
f = sum(f_index);

end

