function p = maxq_10()

p.objective = @maxq;
p.x0 = [1, 1];

end



function f = maxq(x)

n = length(x);
f_index = zeros(n, 1);
for i = 1:n
    f_index(i) = x(i)^2;
end
f = max(f_index);

end