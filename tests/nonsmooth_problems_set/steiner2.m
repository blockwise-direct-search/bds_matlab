function p = steiner2()

p.objective = @steiner2_sub;
A = [
0, 2;
2, 3;
3, -1;
4, -0.5;
5, 2;
6, 2;    
];
a_bar = [5.5, -1];
m = 6;
x0 = ones(2*m, 1);
x0(1) = (A(1, 1) + A(2, 1)) / 3;
for j = 2:m-1
    x0(j) = (x0(j-1) + A(j, 1) + A(j+1, 1)) / 3;
end
x0(m) = (x0(m - 1) + A(m, 1) + a_bar(1)) / 3;
x0(1+m) = (A(1, 2) + A(2, 2)) / 3;
for k = 2:m-1
    x0(k+m) = (x0(k-1+m) + A(k, 2) + A(k+1, 2)) / 3;
end
x0(2*m) = (x0(m-1+m) + A(m, 2) + a_bar(2)) / 3;

p.x0 = x0;

end



function f = steiner2_sub(x)

p = [2, 1, 1, 5, 1, 1];
p_tilde = [1, 1, 2, 3, 2];
A = [
0, 2;
2, 3;
3, -1;
4, -0.5;
5, 2;
6, 2;    
];
a_bar = [5.5, -1];
f_index = zeros(4,1);
m = 6;

for i = 1:length(f_index)
    if i == 1
        f_index(i) = sqrt(x(1)^2 + x(1+m)^2);
    end
    if i == 2
        f_index(i) = sqrt((a_bar(1) - x(m))^2 + (a_bar(2) - x(2*m))^2);
    end

    if i == 3
        for j = 1:m
            f_index(i) = f_index(i) + p(j) * sqrt((x(j) - A(j,1))^2 + (x(j+m) - A(j,2))^2);
        end
    end

    if i == 4
        for j = 1:m-1
            f_index(i) = f_index(i) + p_tilde(j) * sqrt((x(j) - x(j+1))^2 + (x(j+m) - x(j+m+1))^2);
        end
    end

end


f = max(f_index);

end