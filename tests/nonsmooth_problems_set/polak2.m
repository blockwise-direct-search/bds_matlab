function p = polak2()

p.objective = @polak2_sub;
x0 = 0.1*ones(10, 1);
x0(1) = 100;
p.x0 = x0;

end



function f = polak2_sub(x)

f_index = zeros(2,1);

for i = 1:length(f_index)
    if i == 1
        x_tmp = x + eye(10)(:, 2);
        f_index(i) = exp(1e-8*x(1)^2 + x(2:end)'*x(2:end));
    end

    if i == 2
        x_tmp = x - eye(10)(:, 2);
        f_index(i) = exp(1e-8*x(1)^2 + x(2:end)'*x(2:end));
    end
end

f = max(f_index);

end