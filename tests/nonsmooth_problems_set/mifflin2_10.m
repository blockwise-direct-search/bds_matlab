function p = mifflin2_10()

p.objective = @mifflin2;
p.x0 = [1, 1];

end



function f = mifflin2(x)

n = length(x);
f = 0;
for i = 1:n-1
    f = f + (-x(i) + 2*(x(i)^2 + x(i+1)^2 - 1) + 1.75*abs(x(i)^2 + x(i+1)^2 - 1));
end

end