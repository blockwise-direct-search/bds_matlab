function [nb_tag] = get_nb_tag(nb_generator)

solvers_num = length(nb_generator);
nb_tag = strings(1, solvers_num);
for i =1:solvers_num
    if nb_generator(i) == 1
        nb_tag(i) = "one";
    elseif nb_generator(i) == 0.5
        nb_tag(i) = "n";
    elseif nb_generator(i) == 0.25
        nb_tag(i) = "half_n";
    elseif nb_generator(i) == 0.125
        nb_tag(i) = "quarter_n";
    elseif nb_generator(i) == 0.0625
        nb_tag(i) = "half_quarter_n";
    elseif nb_generator(i) == 0
        nb_tag(i) = "none";
    end
end

end

