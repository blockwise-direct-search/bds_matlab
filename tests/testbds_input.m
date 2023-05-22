function [] = testbds_input(parameters)

% Set parameters to an empty structure if it is not supplied.
if nargin < 1
    parameters = struct();
end

%restoredefaultpath;
if isfield(parameters, "matcutest_github_actions")
    if exist(parameters.matcutest_github_actions, "dir")
        addpath(parameters.matcutest_github_actions);
    end
end

% The code of the following lines is for using matcutest.
path_matcutest_server = '/home/htl/local/matcutest/mtools/src';
path_matcutest_local =  '/home/lhtian97/local/matcutest/mtools/src';

if exist(path_matcutest_local, "dir")
    addpath(path_matcutest_local);
end

if exist(path_matcutest_server, "dir")
    addpath(path_matcutest_server);
end

%addpath('/home/lhtian97/bds_new_framework/tests/competitors/prima/matlab/interfaces/');

fullpath = mfilename('fullpath');
[path_tests,~] = fileparts(fullpath);

addpath(path_tests);
cd(path_tests)
% If testdata does not exist, make a new one.
path_testdata = strcat(path_tests, "/testdata");
if ~exist(path_testdata, "dir")
    mkdir(path_testdata);
end

cd ..
path_bds = pwd;
addpath(path_bds);
path_src = strcat(path_bds, "/src");
addpath(path_src);
path_competitors = strcat(path_tests, "/competitors");
addpath(path_competitors);

% The number of solvers. In case that no parameters are input.
if ~isfield(parameters, "num_solvers")
    num_solvers = get_default_testparameters("num_solvers");
else
    num_solvers = parameters.num_solvers;
end

if ~isfield(parameters, "solvers_invoke")
    parameters.solvers_invoke = [];
    for i = 1:num_solvers
        parameters.solvers_invoke = [parameters.solvers_invoke get_default_testparameters("solvers_invoke")];
    end
end

% The code of the following three lines is for running prima first time. If
% we need to compare with prima, then we compile it.
cd(path_tests)
cd ./competitors/prima
setup

if ~isfield(parameters, "solvers_label")
    parameters.solvers_label = [];
    for i = 1:num_solvers
        parameters.solvers_label = [parameters.solvers_label get_default_testparameters("solvers_label")];
    end
end

if ~isfield(parameters, "memory")
    parameters.memory = [];
    for i = 1:num_solvers
        parameters.memory = [parameters.memory get_default_testparameters("memory")];
    end
end

if ~isfield(parameters, "polling_outer")
    parameters.polling_outer = [];
    for i = 1:num_solvers
        parameters.polling_outer = [parameters.polling_outer get_default_testparameters("polling_outer")];
    end
end

if ~isfield(parameters, "polling_inner")
    parameters.polling_inner = [];
    for i = 1:num_solvers
        parameters.polling_inner = [parameters.polling_inner get_default_testparameters("polling_inner")];
    end
end

if ~isfield(parameters, "cycling_inner")
    parameters.cycling_inner = [];
    for i = 1:num_solvers
        parameters.cycling_inner = [parameters.cycling_inner get_default_testparameters("cycling_inner")];
    end
end

if ~isfield(parameters, "solvers_tag")
    parameters.solvers_tag = [];
    for i = 1:num_solvers
        parameters.solvers_tag = [parameters.solvers_tag get_default_testparameters("solvers_tag")];
    end
end

if ~isfield(parameters, "nb_generator")
    parameters.nb_generator = [];
    for i = 1:num_solvers
        parameters.nb_generator = [parameters.nb_generator get_default_testparameters("nb_generator")];
    end
end

if ~isfield(parameters, "nb_tag")
    parameters.nb_tag = [];
    for i = 1:num_solvers
        parameters.nb_tag = [parameters.nb_tag get_default_testparameters("nb_tag")];
    end
end

if ~isfield(parameters, "problems_type")
    parameters.problems_type = get_default_testparameters("problems_type");
end

if ~isfield(parameters, "problems_mindim")
    parameters.problems_mindim = get_default_testparameters("problems_mindim");
end

if ~isfield(parameters, "problems_maxdim")
    parameters.problems_maxdim = get_default_testparameters("problems_maxdim");
end

if ~isfield(parameters, "maxfun_dim")
    parameters.maxfun_dim = get_default_testparameters("maxfun_dim");
    if ~isfield(parameters, "maxfun")
        parameters.maxfun = parameters.maxfun_dim*parameters.problems_maxdim;
    end
end

if ~isfield(parameters, "tau_minimum")
    parameters.tau = 10.^(-1:-1:get_default_testparameters("tau_minimum"));
else
    parameters.tau = 10.^(-1:-1:(-parameters.tau_minimum));
end

if ~isfield(parameters, "direction")
    parameters.direction = [];
    for i = 1:num_solvers
        parameters.direction = [parameters.direction get_default_testparameters("direction")];
    end
end
if isfield(parameters, "parallel")
    parameters.parallel = false;
end
pdfname = "";
% Name pdf automatically (not manually).
for i = 1:num_solvers
    if i > 1
       pdfname = strcat(pdfname, "_"); 
    end
    pdfname = strcat(pdfname, parameters.solvers_tag(i), "_", parameters.nb_tag(i),...
         "_",  num2str(parameters.cycling_inner(i)));
end
pdfname = strcat(pdfname, "_", num2str(parameters.problems_mindim), "_", num2str(parameters.problems_maxdim));
parameters.pdfname = pdfname;
testbds(parameters);
% Delete the path to recover
rmpath(path_tests);
rmpath(path_bds);
rmpath(path_src);
rmpath(path_competitors);

end

