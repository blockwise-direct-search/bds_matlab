function [] = precondition_divide_searching_set(m, nb, cycling, polling)


assert(isintegerscalar(m) && m>0);

assert(isintegerscalar(nb) && nb>0);

assert(isintegerscalar(cycling) && cycling >= 0 && cycling <= 5);

assert(ischstr(polling));


end

