function [nb_tag] = get_nb_tag(nb_generator)

switch nb_generator
    case {1}
        nb_tag = "one";
    case {0.5}
        nb_tag = "n";
    case {0.25}
        nb_tag = "half_n";
    case {0.125}
        nb_tag = "quarter_n";
    case {0.0625}
        nb_tag = "half_quarter_n";
end

end

