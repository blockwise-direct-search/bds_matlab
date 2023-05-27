function [nb_generator] = get_nb(nb_generator)

solvers_num = length(nb_generator);
for i =1:solvers_num
     if ~isnumeric(nb_generator(i))
         if nb_generator(i) == "one"
             nb_generator(i) = 1;
         elseif nb_generator(i) == "n"
             nb_generator(i) = 0.5;
         elseif nb_generator(i) == "half_n"
             nb_generator(i) = 0.25;
         elseif nb_generator(i) == "quarter_n"
             nb_generator(i) = 0.125;
         elseif nb_generator(i) == "half_quarter_n"
             nb_generator(i) = 0.0625;
         end
     end
end

nb_generator = str2double(nb_generator);

end

