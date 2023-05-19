function [output] = testbds(parameters)
% output: performance profile;
% parameters: name; solvers_label; memory; cycling; polling_inner; polling_outer; nb_generator;
% maxfun_dim; maxfun; problems_type; problems_mindim; problems_maxidim; tau;
% clear
% clc

% Get list of problems
s.type = parameters.problems_type; % Unconstrained: 'u'
s.mindim = parameters.problems_mindim; % Minimum of dimension
s.maxdim = parameters.problems_maxdim; % Maximum of dimension
s.blacklist = [];
% Problems that crash.
s.blacklist = [s.blacklist, {}];
% Problems that takes too long to solve.
% {'FBRAIN3LS'} and {'STRATEC'} take too long for fminunc(not for ds and bds).
% {'LRCOVTYPE'} and {'LRIJCNN1'} take long for ds and bds(not for fminunc).
% TODO: check why {'LRIJCNN1'} takes so long to run?
% TODO: check why {'PARKCH'} takes so long to run?
% TODO: check why {'STRATEC'} takes so long to run?
% {'PALMER1C'},{'PALMER2C'},{'PALMER3C'},{'PALMER4C'},{'PALMER5C'},{'PALMER6C'},{'PALMER7C'},{'PALMER8C'}
s.blacklist = [s.blacklist,{'LRCOVTYPE'},{'LRIJCNN1'},{'PARKCH'},{'STRATEC'}...
    ];
problem_names = secup(s);

% list = list(1:min(20,length(list)));
fprintf("We will load %d problems\n\n", length(problem_names))

% Some fixed (relatively) options
% TODO: check what is eps. Read two papers:What Every Computer Scientist Should Know About
% Floating-Point Arithmetic; stability and accuracy numerical(written by Higham).

% Maxfun and maxfun_dim
options_solvers.maxfun = parameters.maxfun; % Maximum of function evaluation
if isfield(parameters, "maxfun_dim")
    options_solvers.maxfun_dim = parameters.maxfun_dim;
end
maxfun = options_solvers.maxfun;

% Termination of stepsize
options_solvers.tol = eps; 
options_solvers.sufficient_decrease_factor = 1e-3;
options_solvers.expand = 2;
options_solvers.shrink = 0.5;

options_solvers.ftarget = -inf;
% TODO: give a reasonable value of the options below.
% TODO: Use consistent way of "" or ''.
options.feature_and_time = "nb";

% acquire fmin and frec
% The difference between solvers_label and name is that solvers_label must be
% different from each other.
options_solvers.solvers = parameters.solvers_invoke;
num_solvers = length(options_solvers.solvers);
num_problems = length(problem_names); % Number of problems
num_random = 1; % Number of random tests(If num_random = 1, it means no random test.)
% The matrix that passed into perfprof.m
frec = NaN(num_problems,num_solvers,num_random,maxfun);
% Store minimum value of the problems of the random test
random_fmin = NaN(num_problems,num_random);

% Some temporary options for test
% noise
options_test.is_noisy = false;
options_test.noise_level = 1e-3;
% relative: (1+noise_level*noise)*f; absolute: f+noise_level*noise
options_test.noise_abs = "relative";
options_test.noise_type = "uniform";

options_test.scaling_matrix = false;
options_test.scaling_matrix_factor = 5;

options_solvers.cycling_inner = parameters.cycling_inner;
options_solvers.polling_inner = parameters.polling_inner;
options_solvers.solvers_label = parameters.solvers_label;
options_solvers.nb_generator = parameters.nb_generator;
options_solvers.memory = parameters.memory;
options_solvers.direction = parameters.direction;

% If parallel is true, use parfor to calculate (parallel computation). ...
% Otherwise, use for to calculate (sequential computation).
if parameters.parallel == true
    parfor i = 1:num_problems
        fval_tmp = NaN(1, num_solvers);
        p = macup(problem_names(1, i));
        for r = 1:num_random
            fprintf('%d(%d). %s\n', i, r, p.name);
            for j = 1:num_solvers
                [fhist,fval] = get_fhist(p, maxfun, j, r, options_solvers, options_test);
                fval_tmp(j) = fval;
                frec(i,j,r,:) = fhist;
            end
            [random_fmin(i,r), I] = min(fval_tmp);
            fprintf('%d %s\n', I, p.name);
        end
    end
else
    for i = 1:num_problems
        fval_tmp = NaN(1, num_solvers);
        p = macup(problem_names(1, i));
        for r = 1:num_random
            fprintf('%d(%d). %s\n', i, r, p.name);
            for j = 1:num_solvers
                [fhist,fval] = get_fhist(p, maxfun, j, r, options_solvers, options_test);
                fval_tmp(j) = fval;
                frec(i,j,r,:) = fhist;
            end
            [random_fmin(i,r), I] = min(fval_tmp);
            fprintf('%d %s\n', I, p.name);
        end
    end
end

if num_random == 1
    fmin = random_fmin;
else
    fmin = min(random_fmin.');
end

% Use time to distinguish
% In matlab, one way to add path is to modify startup.m (Run "which startup.m" to find the location)
time = datestr(now,31);
time = trim_time(time); % trim the form of time
tst = sprintf("test_%s",time); % time stamp to rename
tst = strcat(tst, "_", parameters.pdfname); % rename as mixture of time stamp and parameters
fullpath = mfilename('fullpath');
[path,~] = fileparts(fullpath);
options.path = path;
path_outdir = strcat(options.path, '/testdata/', tst);
options.outdir = strcat(path_outdir, '/perf');
path = strcat(options.path, '/testdata');
% mkdir testdata: test for reproduce this experiment.
cd(path)
mkdir(tst)
% TODO: parellel the code initially (in case experiement is broken)
cd(path_outdir)
mkdir perf
mkdir tests
cd ../../..
path_outdir_src = strcat(path_outdir, '/src');
copyfile('src', sprintf("%s",path_outdir_src));
cd ./tests
path_outdir_tests = strcat(path_outdir, '/tests');
path_outdir_tests_competitors = strcat(path_outdir_tests, '/competitors');
path_outdir_tests_private = strcat(path_outdir_tests, '/private');
copyfile('competitors', sprintf("%s",path_outdir_tests_competitors));
copyfile('private', sprintf("%s",path_outdir_tests_private));
copyfile('testbds.m', sprintf("%s",path_outdir_tests));
copyfile('testbds_parameters.m', sprintf("%s",path_outdir_tests));

% performance profile
tau = parameters.tau; % Tolerance of convergence test in performance profile
tau_length = length(tau);
options_perf.outdir = options.outdir;
options_perf.stamp = time;
options_perf.solvers = parameters.solvers_label;
for l = 1:tau_length
    options_perf.tau = tau(l);
    output = perfprof(frec, fmin, options_perf);
end

cd(options.outdir)
system("pdfunite *.pdf all.pdf");
movefile("all.pdf", sprintf("%s.pdf", parameters.pdfname));
end

