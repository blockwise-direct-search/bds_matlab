function p = maxl()

p.objective = @maxl_sub;
x0 = 1:20;
for i = 1:10
    x0(i+10) = -i;
end

p.x0 = x0;

end



function f = maxl_sub(x)

f = max(abs(x));

end