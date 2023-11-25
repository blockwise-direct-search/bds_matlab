function verify_bds(parameters)
%This function tests the modernized version of Powell's solver against Powell's version, verifying
% whether they produce the same results on CUTEst problems.
%
% 
%
% where
% - `solver` is the name of the solver to test
% - `dimrange` is the vector [mindim, maxdim], or "small", or "big", or "large"
% - `problem` is the name of the problem to test
% - `ir` is the index of the random run in `isequiv`.
% - `nocompile_flag` is either 'nocompile' or 'ncp', indicating not to compile the solves
% - `sequential_flag` (optional) is either 'sequential' or 'seq', which means to test the problems sequentially
% - `reverse_flag` (optional) is either 'reverse' or 'rev', which means to test the solvers in the reverse order
% - `problem_type` can be any of {'u', 'b', 'l', 'n', 'ub', 'ubl', 'ubln', 'bl', 'bln', 'ln'},
%   indicating the problem type to test
%
% Coded by LI Haitian (hai-tian.li@connect.polyu.hk) and Zaikun ZHANG (www.zhangzk.net).
%
% Started: Nov 2023
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Record the current path.
oldpath = path();
% Restore the "right out of the box" path of MATLAB.
restoredefaultpath;
% Record the current directory.
old_dir = pwd();

exception = [];

% Set options to an empty structure if it is not provided.
if nargin < 1
    parameters = struct();
end

try

    % Compile the version of norma.
    path_norma = locate_norma();
    path_verify_bds = fileparts(mfilename('fullpath'));
    cd(path_norma{1});
    setup
    cd(path_verify_bds);

    % Compile the version of modern repository.
    path_root = fileparts(path_verify_bds);
    cd(path_root);
    setup
    cd(path_verify_bds);

    solvers = {"bds", "bds_norma"};

    % Tell MATLAB where to find MatCUTEst.
    locate_matcutest();

    % Record `olddir` in `options` so that we can come back to `olddir` during `isequiv` if
    % necessary (for example, when a single test fails).
    parameters.olddir = old_dir;

    % Get list of problems
    if isfield(parameters, "problem_type")
        s.type = parameters.problem_type;
    else
        s.type = get_default_profile_options("problem_type");
    end

    if isfield(parameters, "problem_mindim")
        s.mindim = parameters.problem_mindim;
    else
        s.mindim = get_default_profile_options("problem_mindim");
    end

    if isfield(parameters, "problem_maxdim")
        s.maxdim = parameters.problem_maxdim;
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

    if ~isfield(parameters, "problem_names")
        problem_names = secup(s);
    else
        problem_names = parameters.problem_names;
    end
    
    fprintf("We will load %d problems\n\n", length(problem_names))

    % Get the number of problems.
    num_problems = length(problem_names);

    if isfield(parameters, 'prec')
        prec = parameters.prec;
    else
        prec = 0;
    end

    if ~isfield(parameters, "parallel")
        parallel = get_default_profile_options("parallel");
    else
        parallel = parameters.parallel;
    end
    parameters.sequential = ~parallel;

    if ~isfield(parameters, "i_run_init")
        i_run_init = 1;
    else
        i_run_init = parameters.i_run_init;
    end

    if ~isfield(parameters, "num_random")
        num_random = 20;
    else
        num_random = parameters.num_random;
    end

    if isfield(parameters, 'single_test')
        single_test = parameters.single_test;
    else
        single_test = (num_random == 1);
    end
    
    if parallel
        parfor i_problem = 1:num_problems
            p = macup(problem_names(1, i_problem));
            for i_run = i_run_init:i_run_init+num_random-1
                fprintf("%d(%d). %s\n", i_problem, i_run, p.name);
                iseqiv(solvers, p, i_run, single_test, prec, parameters);
            end
        end
    else
        for i_problem = 1:num_problems
            p = macup(problem_names(1, i_problem));
            for i_run = i_run_init:i_run_init+num_random-1
                fprintf("%d(%d). %s\n", i_problem, i_run, p.name);
                iseqiv(solvers, p, i_run, single_test, prec, parameters);
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