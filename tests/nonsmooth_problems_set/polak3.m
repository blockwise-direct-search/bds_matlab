function p = polak3()

p.objective = @polak3_sub;
p.x0 = ones(11, 1);

end



function f = polak3_sub(x)

f_index = zeros(10,1);

for i = 1:length(f_index)
    for j = 0:10
        f_index(i) = f_index(i) + (1/(i+j))*exp((x(j+1) - sin(i - 1 + 2*j))^2);
    end
end

f = max(f_index);

end