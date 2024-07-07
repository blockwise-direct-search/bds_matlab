function p = gill()

p.objective = @gill_sub;
p.x0 = -0.1*zeros(10, 1);

end



function f = gill_sub(x)

f_index = zeros(3, 1);
for i = 1:10
    f_index(1) = f_index(1) + (x(i) - 1)^2 + 1e-3*(x(i)^2 - 0.25)^2;
end

for i = 2:30
    tmp1 = 0;
    for j = 2:10
        tmp1 = tmp1 + x(j)*(j - 1) * ((i - 1)/29)^(j - 2);
    end

    tmp2 = 0;
    for j = 1:10
        tmp2 = tmp2 + x(j) * ((i - 1)/29)^(j - 1); 
    end
    tmp2 = tmp2^2;

    f_index(2) = f_index(2) + (tmp1 - tmp2 - 1)^2 + x(1)^2 + (x(2) - x(1)^2 - 1)^2;
end


for i = 2:10
    f_index(3) = f_index(3) + 100*(x(i) - x(i-1)^2)^2 + (x(i-1) - 1)^2;
end

f = max(f_index);

end