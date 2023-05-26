function [options] = get_options(p, j, name_solver, options_solvers, options)

prima_list = ["cobyla", "uobyqa", "newuoa", "bobyqa", "lincoa"];
if ~isempty(find(prima_list == name_solver, 1))
    name_solver = "prima";
end

maxfun = options.maxfun;

if name_solver == "blockwise_direct_search"
    
    % Polling strategies should be defined in the loop!!!
    options.polling_inner = options_solvers.polling_inner(j);
    
    % Strategy of blocking
    % If nb_generator<1, nb may be flexible by different
    % dimensions, otherwise nb is fixed.
    % 2.5 is warning!
    x0 = p.x0;
    dim = length(x0);
    if options_solvers.nb_generator(j) >= 1
        if ceil(options_solvers.nb_generator(j)) == options_solvers.nb_generator(j)
            options.nb = options_solvers.nb_generator(j);
        else
            options.nb = ceil(options_solvers.nb_generator(j));
            disp("Wrong input of nb_generator");
        end
    else
        options.nb = ceil(2*dim*options_solvers.nb_generator(j));
    end
    
    % Strategy of memory, cycling and polling_inner (Memory vs Nonmemory when cycling)
    options.memory = options_solvers.memory(j);
    options.cycling_inner = options_solvers.cycling_inner(j);
    options.direction = options_solvers.direction(j);
    options.blocks_strategy = options_solvers.blocks_strategy(j);
    
elseif name_solver == "bds_polling"
    
    % Polling strategies should be defined in the loop!!!
    options.polling_inner = options_solvers.polling_inner(j);
    
    % Strategy of blocking
    % If nb_generator<1, nb may be flexible by different
    % dimensions, otherwise nb is fixed.
    % 2.5 is warning!
    x0 = p.x0;
    dim = length(x0);
    if options_solvers.nb_generator(j) >= 1
        if ceil(options_solvers.nb_generator(j)) == options_solvers.nb_generator(j)
            options.nb = options_solvers.nb_generator(j);
        else
            options.nb = ceil(options_solvers.nb_generator(j));
            disp("Wrong input of nb_generator");
        end
    else
        options.nb = ceil(2*dim*options_solvers.nb_generator(j));
    end
    
    % Strategy of memory, cycling and polling_inner (Memory vs Nonmemory when cycling)
    options.memory = options_solvers.memory(j);
    options.cycling_inner = options_solvers.cycling_inner(j);
    options.direction = options_solvers.direction(j);    
    
elseif name_solver == "prima"
    options.output_xhist = true;
    % An indicator: it can attain 0, 1, 2, 3, -1, -2, -3. Default value is
    % 0. More absolute value of iprint, more information will be printed on command
    % window. When the value of iprint is negative, no information will be
    % printed on command window and will be stored in a file.
    options.iprint = 0;
    % options.classical = true;
    
elseif name_solver == "matlab_fminsearch"
    oldopts = optimset('MaxFunEvals', maxfun);
    options = optimset(oldopts, 'MaxIter', maxfun);
    
elseif name_solver == "matlab_fminunc"    
    options = optimoptions(@fminunc,'MaxFunctionEvaluations', maxfun, 'MaxIterations', maxfun, 'StepTolerance', 1e-12);
    
else
    disp("there are no options for the j-th solver");
end


end

