function [nb_generator] = get_nb(nb_tag)

switch nb_tag
    case {"one"}
        nb_generator = 1;
    case {"n"}
        nb_generator = 0.5;
    case {"half_n"}
        nb_generator = 0.25;
    case {"quarter_n"}
        nb_generator = 0.125;
    case {"half_quarter_n"}
        nb_generator = 0.0625;
end

end

