function p = watson()

p.objective = @watson_sub;
p.x0 = zeros(20, 1);

end



function f = watson_sub(x)

f_index = zeros(31, 1);
f_index(1) = x(1);
f_index(2) = x(2) - x(1)^2 - 1;
for i = 3:31
    for j = 1:20
        tmp1 = 0;
        for k = 2:20
            tmp1 = tmp1 + (k -1)*x(k)((i - 2)/29)^(k-2);
        end

        tmp2 = 0;
        for l = 1:20
            tmp2 = tmp2 + x(l)*((i - 2)/29)^(l-1);
        end
        tmp2 = tmp2^2;

    end
    f_index(i) = tmp1 - tmp2;
end

f = max(abs(f_index));

end