function [output] = profile_bds(parameters)
% output: performance profile;
% parameters: name; solvers_legend; memory; cycling; polling_inner; polling_outer; nb_generator;
% maxfun_dim; maxfun; problems_type; problems_mindim; problems_maxidim; tau;
% clear
% clc

% Set parameters to an empty structure if it is not supplied.
if nargin < 1
    parameters = struct();
end

restoredefaultpath;

% % The code of the following lines is for using matcutest.
% path_matcutest_server = '/home/htl/local/matcutest/mtools/src';
% path_matcutest_local =  '/home/lhtian97/local/matcutest/mtools/src';
%addpath('/home/lhtian97/bds_new_framework/tests/competitors/prima/matlab/interfaces/');

fullpath = mfilename("fullpath");
[path_tests,~] = fileparts(fullpath);
parameters.path_tests = path_tests;
path_bds = fileparts(path_tests);
parameters.path_bds = path_bds;
addpath(path_tests);
% If testdata does not exist, make a new one.
path_testdata = strcat(path_tests, "/testdata");
if ~exist(path_testdata, "dir")
    mkdir(path_testdata);
end

addpath(path_bds);
path_src = fullfile(path_bds, "src");
parameters.path_src = path_src;
addpath(path_src);
path_competitors = fullfile(path_tests, "competitors");
addpath(path_competitors);
parameters.path_competitors = path_competitors;
path_competitors_mnewuoa = fullfile(path_competitors, "mnewuoa");
addpath(path_competitors_mnewuoa);
path_competitors_matlab_functions = fullfile(path_competitors, "matlab_functions");
addpath(path_competitors_matlab_functions);

assert(isfield(parameters, "solvers_invoke"));

parameters = get_profile_options(parameters);

% Tell MATLAB where to find MatCUTEst.
locate_matcutest();
% Tell MATLAB where to find prima.
locate_prima();

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
% s.blacklist = [s.blacklist,{'LRCOVTYPE'},{'LRIJCNN1'},{'PARKCH'},{'STRATEC'}...
%     ];
s.blacklist = [{'LRCOVTYPE'},{'LRIJCNN1'},{'PARKCH'},{'STRATEC'}];
%    {'GAUSS1LS'}, {'GAUSS2LS'}, {'GAUSS3LS'}, {'HEART8LS'}, ...
%     {'PALMER1C'}, {'PALMER2C'}, {'PALMER3C'}, {'PALMER4C'}, {'PALMER6C'}, ...
%     {'PALMER7C'}, {'PALMER8C'}, {'VESUVIALS'}, {'VESUVIOLS'}, {'VESUVIOULS'}, ...
%     {'VIBRBEAM'}

problem_names = secup(s);

% list = list(1:min(20,length(list)));
fprintf("We will load %d problems\n\n", length(problem_names))

% Some fixed (relatively) options
% TODO: check what is eps. Read two papers:What Every Computer Scientist Should Know About
% Floating-Point Arithmetic; stability and accuracy numerical(written by Higham).

% Maxfun and maxfun_dim
solver_options.maxfun = parameters.maxfun; % Maximum of function evaluation
if isfield(parameters, "maxfun_dim")
    solver_options.maxfun_dim = parameters.maxfun_dim;
end
maxfun = solver_options.maxfun;

% Parameters of stepsize
solver_options.StepTolerance = parameters.StepTolerance;
solver_options.sufficient_decrease_factor = parameters.sufficient_decrease_factor;
solver_options.expand = parameters.expand;
solver_options.shrink = parameters.shrink;
solver_options.alpha_init = parameters.alpha_init;

if isfield(parameters, "powell_factor")
    solver_options.powell_factor = parameters.powell_factor;
end

if isfield(parameters, "accept_simple_decrease")
    solver_options.accept_simple_decrease = parameters.accept_simple_decrease;
end

% Parameters of ftarget
solver_options.ftarget = parameters.ftarget;

% Initialize fmin and frec
% The difference between solvers_legend and name is that solvers_legend must be
% different from each other.
solver_options.solvers = parameters.solvers_invoke;
num_solvers = length(solver_options.solvers);
num_problems = length(problem_names); % Number of problems
num_random = parameters.num_random; % Number of random tests(If num_random = 1, it means no random test.)
% The matrix that passed into perfprof.m
frec = NaN(num_problems,num_solvers,num_random,maxfun);
% Store minimum value of the problems of the random test
fmin = NaN(num_problems, num_random);

% Some temporary options for test
% noise
test_options.is_noisy = parameters.is_noisy;
test_options.noise_level = parameters.noise_level;
% relative: (1+noise_level*noise)*f; absolute: f+noise_level*noise
test_options.is_abs_noise = parameters.is_abs_noise;
test_options.noise_type = parameters.noise_type;
test_options.num_random = parameters.num_random;

test_options.scaling_matrix = false;
test_options.scaling_matrix_factor = 5;

solver_options.fminunc_type = parameters.fminunc_type;

if isfield(parameters, "blocks_strategy")
    solver_options.blocks_strategy = parameters.blocks_strategy;
end

if isfield(parameters, "nb_generator")
    solver_options.nb_generator = parameters.nb_generator;
end

if isfield(parameters, "randomized_strategy")
    solver_options.randomized_strategy = parameters.randomized_strategy;
end

solver_options.cycling_inner = parameters.cycling_inner;
solver_options.polling_inner = parameters.polling_inner;
solver_options.solvers_legend = parameters.solvers_legend;
solver_options.with_memory = parameters.with_memory;
solver_options.direction = parameters.direction;

% If parallel is true, use parfor to calculate (parallel computation). ...
% Otherwise, use for to calculate (sequential computation).
if parameters.parallel == true
    parfor i = 1:num_problems
        fval_tmp = NaN(1, num_solvers);
        p = macup(problem_names(1, i));
        for r = 1:num_random
            fprintf("%d(%d). %s\n", i, r, p.name);
            for j = 1:num_solvers
                [fhist,fval] = get_fhist(p, maxfun, j, r, solver_options, test_options);
                fval_tmp(j) = fval
                frec(i,j,r,:) = fhist;
            end
             [fmin(i,r), I] = min(fval_tmp);
            index_min = find(fval_tmp <= fmin(i,r));
            fprintf("%s %s\n", sprintf('%d ', index_min), p.name);
        end
    end
else
    for i = 1:num_problems
        fval_tmp = NaN(1, num_solvers);
        p = macup(problem_names(1, i));
        for r = 1:num_random
            fprintf("%d(%d). %s\n", i, r, p.name);
            for j = 1:num_solvers
                [fhist,fval] = get_fhist(p, maxfun, j, r, solver_options, test_options);
                fval_tmp(j) = fval;
                frec(i,j,r,:) = fhist;
            end
            [fmin(i,r), I] = min(fval_tmp);
            index_min = find(fval_tmp <= fmin(i,r));
            fprintf("%s %s\n", sprintf('%d ', index_min), p.name);
        end
    end
end

% Use time to distinguish
% In matlab, one way to add path is to modify startup.m (Run "which startup.m" to find the location)
time = datestr(now,31);
time = trim_time(time); % trim the form of time
tst = sprintf("test_%s",time); % time stamp to rename
tst = strcat(tst, "_", parameters.pdfname); % rename as mixture of time stamp and parameters
options.path = parameters.path_tests;
path_testdata = fullfile(options.path, "testdata");
path_testdata_outdir = fullfile(options.path, "testdata", tst);

% mkdir a new folder to save numerical results and source code.
mkdir(path_testdata, tst);
mkdir(path_testdata_outdir, "perf");
options.outdir = fullfile(path_testdata_outdir, "perf");
mkdir(path_testdata_outdir, "src");
path_testdata_src = fullfile(path_testdata_outdir, "src");
mkdir(path_testdata_outdir, "tests");
path_testdata_tests = fullfile(path_testdata_outdir, "tests");
path_testdata_competitors = fullfile(path_testdata_tests, "competitors");
mkdir(path_testdata_competitors);
path_testdata_private = fullfile(path_testdata_tests, "private");
mkdir(path_testdata_private);

% Copy the source code and test code to path_outdir.
copyfile(fullfile(parameters.path_src, "*"), path_testdata_src);
parameters.path_competitors = fullfile(parameters.path_tests, "competitors");
copyfile(fullfile(parameters.path_competitors, "*"), path_testdata_competitors);
copyfile(fullfile(parameters.path_tests, "private", "*"), path_testdata_private);

source_folder = parameters.path_tests;
destination_folder = path_testdata_tests;

% Get all files in the source folder.
file_list = dir(fullfile(source_folder, '*.*'));
file_list = file_list(~[file_list.isdir]);

% Copy all files (excluding subfolders) to the destination folder.
for i = 1:numel(file_list)
    source_file = fullfile(source_folder, file_list(i).name);
    destination_file = fullfile(destination_folder, file_list(i).name);
    copyfile(source_file, destination_file);
end

% performance profile
tau = parameters.tau; % Tolerance of convergence test in performance profile
tau_length = length(tau);
options_perf.outdir = options.outdir;
options_perf.stamp = time;
options_perf.solvers = parameters.solvers_legend;
options_perf.natural_stop = false;
for l = 1:tau_length
    options_perf.tau = tau(l);
    output = perfprof(frec, fmin, options_perf);
end

path_temporary = pwd;
cd(options.outdir);

% Initialize string variable.
pdfFiles = dir(fullfile(options.outdir, '*.pdf'));

% Store filename in a cell.
pdfNamesCell = cell(numel(pdfFiles), 1);
for i = 1:numel(pdfFiles)
    pdfNamesCell{i} = pdfFiles(i).name;
end

% Use the strjoin function to concatenate the elements in a cell array into a single string.
inputfiles = strjoin(pdfNamesCell, ' ');

% Remove spaces at the beginning of a string.
inputfiles = strtrim(inputfiles);

% Merge pdf.
outputfile = 'all.pdf';
system(['bash ', fullfile(parameters.path_tests, 'private', 'compdf'), ' ', inputfiles, ' -o ', outputfile]);
% Rename pdf
movefile("all.pdf", sprintf("%s.pdf", parameters.pdfname));

% Delete the path to recover
rmpath(path_tests);
rmpath(path_bds);
rmpath(path_src);
rmpath(path_competitors);
rmpath(path_competitors_mnewuoa);
rmpath(path_competitors_matlab_functions);

% Restore the path
cd(path_temporary);

end
