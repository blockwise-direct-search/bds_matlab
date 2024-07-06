function p = lq_10()

p.objective = @lq;
p.x0 = [1, 1];

end



function f = lq(x)

n = length(x);
f_index = zeros(n-1, 1);
for i = 1:n-1
    f_index(i) = max(-x(i)-x(i+1), -x(i)-x(i+1)+(x(i)^2 + x(i+1)^2 - 1)); 
end
f = sum(f_index);

end

