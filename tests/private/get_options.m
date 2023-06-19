function [options] = get_options(p, j, name_solver, solver_options, options)

bds_list = ["bds", "bds_powell"];
prima_list = ["cobyla", "uobyqa", "newuoa", "bobyqa", "lincoa"];

maxfun = options.maxfun;

if any(contains(bds_list, name_solver, 'IgnoreCase', true))

    % Polling strategies should be defined in the loop!!!
    options.polling_inner = solver_options.polling_inner(j);

    % Strategy of blocking
    % If nb_generator<1, nb may be flexible by different
    % dimensions, otherwise nb is fixed.
    % 2.5 is warning!
    x0 = p.x0;
    dim = length(x0);
    if solver_options.nb_generator(j) >= 1
        if ceil(solver_options.nb_generator(j)) == solver_options.nb_generator(j)
            options.nb = solver_options.nb_generator(j);
        else
            options.nb = ceil(solver_options.nb_generator(j));
            disp("Wrong input of nb_generator");
        end
    else
        options.nb = ceil(2*dim*solver_options.nb_generator(j));
    end

    % Strategy of with_memory, cycling and polling_inner (Memory vs Nonwith_memory when cycling)
    options.with_memory = solver_options.with_memory(j);
    options.cycling_inner = solver_options.cycling_inner(j);
    options.direction = solver_options.direction(j);

    % Options of step size
    options.StepTolerance = solver_options.StepTolerance;
    options.sufficient_decrease_factor = solver_options.sufficient_decrease_factor;
    options.expand = solver_options.expand;
    options.shrink = solver_options.shrink;
    options.alpha_init = solver_options.alpha_init;

    if isfield(solver_options, "Algorithm")
        options.Algorithm = solver_options.Algorithm(j);
    end
    if isfield(solver_options, "powell_factor")
        options.powell_factor = solver_options.powell_factor(j);
    end

    if isfield(solver_options, "accept_simple_decrease")
        options.accept_simple_decrease = solver_options.accept_simple_decrease(j);
    end

elseif any(contains(prima_list, name_solver, 'IgnoreCase', true))
    options.output_xhist = true;
    % An indicator: it can attain 0, 1, 2, 3, -1, -2, -3. Default value is
    % 0. More absolute value of iprint, more information will be printed on command
    % window. When the value of iprint is negative, no information will be
    % printed on command window and will be stored in a file.
    options.iprint = 0;

    if isfield(solver_options, "classical")
        options.classical = solver_options.classical;
    end

    % Options of trust region radius
    options.rhobeg = solver_options.alpha_init;
    options.rhoend = solver_options.StepTolerance;

elseif name_solver == "matlab_fminsearch"
    options = optimset('MaxFunEvals', maxfun, 'maxiter', maxfun, 'tolfun',...
        solver_options.StepTolerance, 'tolx', solver_options.StepTolerance);

elseif name_solver == "matlab_fminunc"
    options = optimoptions('fminunc', 'Algorithm', 'quasi-newton', ...
        'HessUpdate', solver_options.fminunc_type, 'MaxFunctionEvaluations',...
    maxfun, 'MaxIterations', maxfun, 'ObjectiveLimit', solver_options.ftarget,...
    'StepTolerance', solver_options.StepTolerance, 'OptimalityTolerance', solver_options.StepTolerance);

elseif name_solver == "matlab_patternsearch"
    options = optimoptions('patternsearch','MaxIterations', maxfun,...
    'MaxFunctionEvaluations', maxfun, 'FunctionTolerance', solver_options.StepTolerance,...
        'TolMesh', solver_options.StepTolerance, 'StepTolerance', solver_options.StepTolerance);

else
    fprintf("%s\n", name_solver)
    disp("there are no options for the j-th solver");
end


end
