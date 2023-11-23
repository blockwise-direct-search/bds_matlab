function verify(parameters)

% Record the current path.
oldpath = path();
% Restore the "right out of the box" path of MATLAB.
restoredefaultpath;
% Record the current directory.
old_dir = pwd();

exception = [];

try

    % Compile the version of norma.
    path_norma = locate_norma();
    path_verify = fileparts(mfilename('fullpath'));
    cd(path_norma{1});
    setup
    cd(path_verify);

    % Compile the version of modern repository.
    path_root = fileparts(path_verify);
    cd(path_root);
    setup
    cd(path_verify);

    solvers = {"bds", "bds_norma"};

    % Tell MATLAB where to find MatCUTEst.
    locate_matcutest();

    % Record `olddir` in `options` so that we can come back to `olddir` during `isequiv` if
    % necessary (for example, when a single test fails).
    parameters.olddir = old_dir;

    % Get list of problems
    if isfield(parameters, "problems_type")
        s.type = parameters.problems_type;
    else
        s.type = get_default_profile_options("problems_type");
    end

    if isfield(parameters, "problems_mindim")
        s.mindim = parameters.problems_mindim;
    else
        s.mindim = get_default_profile_options("problems_mindim");
    end

    if isfield(parameters, "problems_maxdim")
        s.maxdim = parameters.problems_maxdim;
    else
        s.maxdim = 20;
    end
    
    s.blacklist = [];

    if s.mindim >= 6
        s.blacklist = [parameters.blacklist, { 'ARGTRIGLS', 'BROWNAL', ...
            'COATING', 'DIAMON2DLS', 'DIAMON3DLS', 'DMN15102LS', ...
            'DMN15103LS', 'DMN15332LS', 'DMN15333LS', 'DMN37142LS', ...
            'DMN37143LS', 'ERRINRSM', 'HYDC20LS', 'LRA9A', ...
            'LRCOVTYPE', 'LUKSAN12LS', 'LUKSAN14LS', 'LUKSAN17LS', 'LUKSAN21LS', ...
            'LUKSAN22LS', 'MANCINO', 'PENALTY2', 'PENALTY3', 'VARDIM',
            }];
    end

    problem_names = secup(s);

    fprintf("We will load %d problems\n\n", length(problem_names))

    % Get the number of problems.
    num_problems = length(problem_names);

    if isfield(parameters, 'prec')
        prec = parameters.prec;
    else
        prec = 0;
    end

    if isfield(parameters, 'single_test')
        single_test = parameters.single_test;
    else
        single_test = false;
    end

    if ~isfield(parameters, "parallel")
        parallel = get_default_profile_options("parallel");
    end

    if ~isfield(parameters, "num_random")
        num_random = 20;
    else
        num_random = parameters.num_random;
    end
    
    if parallel
        parfor i_problem = 1:num_problems
            p = macup(problem_names(1, i_problem));
            for i_run = 1:num_random
                fprintf("%d(%d). %s\n", i_problem, i_run, p.name);
                iseqiv(solvers, p, i_run, prec, single_test, parameters);
            end
        end
    else
        for i_problem = 1:num_problems
            p = macup(problem_names(1, i_problem));
            for i_run = 1:num_random
                fprintf("%d(%d). %s\n", i_problem, i_run, p.name);
                iseqiv(solvers, p, i_run, prec, single_test, parameters);
            end
        end
    end

catch exception

    % Do nothing for the moment.

end

% Restore the path to oldpath.
setpath(oldpath);
% Go back to the original directory.
cd(old_dir);

if ~isempty(exception)  % Rethrow any exception caught above.
    rethrow(exception);
end

end