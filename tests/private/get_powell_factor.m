function powell_factor = get_powell_factor(powell_factor_level)

switch powell_factor_level
    case {"one"}
        powell_factor = 0.1;
    case {"two"}
        powell_factor = 1e-2;
    case {"three"}
        powell_factor = 1e-3;
    case {"four"}
        powell_factor = 1e-4;
    case {"five"}
        powell_factor = 1e-5;
    case {"six"}
        powell_factor = 1e-6;
    case {"seven"}
        powell_factor = 1e-7;
    case {"eight"}
        powell_factor = 1e-8;
    case {"nine"}
        powell_factor = 1e-9;
    case {"ten"}
        powell_factor = 1e-10;
end

end

