function [output] = profile_bds(parameters)
% Draw performance profiles.

% Record the current path.
oldpath = path(); 
% Restore the "right out of the box" path of MATLAB. 
restoredefaultpath;  
% Record the current directory.
old_dir = pwd();

% Add the paths that we need to use in the performance profile into the MATLAB
% search path.
fullpath = mfilename("fullpath");
[path_tests,~] = fileparts(fullpath);
parameters.path_tests = path_tests;
path_bds = fileparts(path_tests);
parameters.path_bds = path_bds;
addpath(path_tests);

% If the folder of testdata does not exist, make a new one.
path_testdata = fullfile(path_tests, "testdata");
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

% In case no solvers are input, then throw an error.
if ~isfield(parameters, "solvers_invoke")
    error("There should be at least two solvers."); 
end

% Get the parameters that the test needs.
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
% TODO: check what is eps. Read two papers: What Every Computer Scientist Should Know About
% Floating-Point Arithmetic; stability and accuracy numerical(written by Higham).

% Set maxfun and maxfun_dim.
solver_options.maxfun = parameters.maxfun;
if isfield(parameters, "maxfun_dim")
    solver_options.maxfun_dim = parameters.maxfun_dim;
end
maxfun = solver_options.maxfun;

% Set parameters of stepsize.
solver_options.StepTolerance = parameters.StepTolerance;
solver_options.sufficient_decrease_factor = parameters.sufficient_decrease_factor;
solver_options.forcing_function = parameters.forcing_function;
solver_options.expand = parameters.expand;
solver_options.shrink = parameters.shrink;
solver_options.alpha_init = parameters.alpha_init;
solver_options.num_random_vectors = parameters.num_random_vectors;

if isfield(parameters, "powell_factor")
    solver_options.powell_factor = parameters.powell_factor;
end

if isfield(parameters, "accept_simple_decrease")
    solver_options.accept_simple_decrease = parameters.accept_simple_decrease;
end

if isfield(parameters, "shuffling_period")
    solver_options.shuffling_period = parameters.shuffling_period;
end

if isfield(parameters, "replacement_delay")
    solver_options.replacement_delay = parameters.replacement_delay;
end

% Set ftarget of objective function.
solver_options.ftarget = parameters.ftarget;

% Initialize fmin and frec.
solver_options.solvers = parameters.solvers_invoke;
num_solvers = length(solver_options.solvers);
% Get number of problems.
num_problems = length(problem_names); 
% Get Number of random tests(If num_random = 1, it means no random test).
num_random = parameters.num_random; 
% Record minimum value of the problems of the random test.
if parameters.is_noisy && strcmpi(parameters.fmin_type, "real-randomized")
    fmin = NaN(num_problems, num_random+1);
    % Get the matrix that passed into perfprof.m. In this case, test results without noise
    % also need to be passed.
    frec = NaN(num_problems, num_solvers, num_random+1, maxfun);
else
    % Get the matrix that passed into perfprof.m.
    fmin = NaN(num_problems, num_random);
    frec = NaN(num_problems, num_solvers, num_random, maxfun);
end 

% Set noisy parts of test.
test_options.is_noisy = parameters.is_noisy;
test_options.noise_level = parameters.noise_level;
% Relative: (1+noise_level*noise)*f; absolute: f+noise_level*noise
test_options.is_abs_noise = parameters.is_abs_noise;
test_options.noise_type = parameters.noise_type;
test_options.num_random = parameters.num_random;

% Set scaling matrix.
test_options.scaling_matrix = false;
test_options.scaling_matrix_factor = 5;

solver_options.fminunc_type = parameters.fminunc_type;

if isfield(parameters, "Algorithm")
    solver_options.Algorithm = parameters.Algorithm;
end

if isfield(parameters, "nb_generator")
    solver_options.nb_generator = parameters.nb_generator;
end

if isfield(parameters, "randomized_strategy")
    solver_options.randomized_strategy = parameters.randomized_strategy;
end

if isfield(parameters, "classical")
    solver_options.classical = parameters.classical;
end

if isfield(parameters, "cunxin_factor")
    solver_options.cunxin_factor = parameters.cunxin_factor;
end

if isfield(parameters, "cunxin_factor_period")
    solver_options.cunxin_factor_period = parameters.cunxin_factor_period;
end

solver_options.cycling_inner = parameters.cycling_inner;
solver_options.polling_inner = parameters.polling_inner;
solver_options.solvers_legend = parameters.solvers_legend;
solver_options.with_cycling_memory = parameters.with_cycling_memory;
solver_options.direction = parameters.direction;


% If parameters.noise_initial_point is true, then initial point will be 
% selected for each problem num_random times.
% parameters.fmintype is set to be "randomized" defaultly, then there is
% no need to test without noise, which makes the curve of performance profile
% more higher. If parallel is true, use parfor to calculate (parallel computation), 
% otherwise, use for to calculate (sequential computation).
if parameters.parallel == true
    parfor i = 1:num_problems
        fval_tmp = NaN(1, num_solvers);
        p = macup(problem_names(1, i));
        for r = 1:num_random
            if parameters.noise_initial_point
                dim = length(p.x0);
                rr = randn(dim, 1);
                rr = rr / norm(rr);
                p.x0 = p.x0 + 10 * max(1, norm(p.x0)) * rr;
            end
            fprintf("%d(%d). %s\n", i, r, p.name);
            for j = 1:num_solvers
                [fhist,fval] = get_fhist(p, maxfun, j, r, solver_options, test_options);
                fval_tmp(j) = fval
                frec(i,j,r,:) = fhist;
            end
            [fmin(i,r), ~] = min(fval_tmp);
            index_min = find(fval_tmp <= fmin(i,r));
            fprintf("%s %s\n", sprintf('%d ', index_min), p.name);
        end
    end
else
    for i = 1:num_problems
        fval_tmp = NaN(1, num_solvers);
        p = macup(problem_names(1, i));
        for r = 1:num_random
            if parameters.noise_initial_point
                dim = length(p.x0);
                rr = randn(dim, 1);
                rr = rr / norm(rr);
                p.x0 = p.x0 + 10 * max(1, norm(p.x0)) * rr;
            end
            fprintf("%d(%d). %s\n", i, r, p.name);
            for j = 1:num_solvers
                [fhist,fval] = get_fhist(p, maxfun, j, r, solver_options, test_options);
                fval_tmp(j) = fval;
                frec(i,j,r,:) = fhist;
            end
            [fmin(i,r), ~] = min(fval_tmp);
            index_min = find(fval_tmp <= fmin(i,r));
            fprintf("%s %s\n", sprintf('%d ', index_min), p.name);
        end
    end
end

% If parameters.fmintype = "real-randomized", then test without noise
% should be conducted and fmin might be smaller, which makes curves
%  of performance profile more lower.
if parameters.is_noisy && strcmpi(parameters.fmin_type, "real-randomized")
    test_options.is_noisy = false;
    r = num_random+1;
    if parameters.parallel == true
        parfor i = 1:num_problems
            fval_tmp = NaN(1, num_solvers);
            p = macup(problem_names(1, i));
            if parameters.noise_initial_point
                dim = length(p.x0);
                rr = randn(dim, 1);
                rr = rr / norm(rr);
                p.x0 = p.x0 + 10 * max(1, norm(p.x0)) * rr;
            end
            fprintf("%d(%d). %s\n", i, r, p.name);
            for j = 1:num_solvers
                [fhist,fval] = get_fhist(p, maxfun, j, r, solver_options, test_options);
                fval_tmp(j) = fval;
                frec(i,j,r,:) = fhist;
            end
            [fmin(i,r), ~] = min(fval_tmp);
            index_min = find(fval_tmp <= fmin(i,r));
            fprintf("%s %s\n", sprintf('%d ', index_min), p.name);
        end
    else
        for i = 1:num_problems
            fval_tmp = NaN(1, num_solvers);
            p = macup(problem_names(1, i));
            if parameters.noise_initial_point
                dim = length(p.x0);
                rr = randn(dim, 1);
                rr = rr / norm(rr);
                p.x0 = p.x0 + 10 * max(1, norm(p.x0)) * rr;
            end
            fprintf("%d(%d). %s\n", i, r, p.name);
            for j = 1:num_solvers
                [fhist,fval] = get_fhist(p, maxfun, j, r, solver_options, test_options);
                fval_tmp(j) = fval;
                frec(i,j,r,:) = fhist;
            end
            [fmin(i,r), ~] = min(fval_tmp);
            index_min = find(fval_tmp <= fmin(i,r));
            fprintf("%s %s\n", sprintf('%d ', index_min), p.name);
        end    
    end
end 

% Use time to distinguish.
time = datetime("now");
time_str = sprintf('%04d-%02d-%02d %02d:%02d:%02.0f', year(time), ...
    month(time), day(time), hour(time), minute(time), second(time));
% Trim time string.
time_str = trim_time(time_str); 
tst = sprintf("test_%s", time_str); 
% Rename tst as mixture of time stamp and pdfname.
tst = strcat(tst, "_", parameters.pdfname); 
options.path = parameters.path_tests;
path_testdata = fullfile(options.path, "testdata");
path_testdata_outdir = fullfile(options.path, "testdata", tst);

% Make a new folder to save numerical results and source code.
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

% Draw performance profiles.
% Set tolerance of convergence test in performance profile.
tau = parameters.tau; 
tau_length = length(tau);
options_perf.outdir = options.outdir;
options_perf.stamp = time_str;
options_perf.solvers = parameters.solvers_legend;
options_perf.natural_stop = false;
if parameters.is_noisy && strcmpi(parameters.fmin_type, "real-randomized")
    fmin = min(fmin');
end  
for l = 1:tau_length
    options_perf.tau = tau(l);
    output = perfprof(frec, fmin, options_perf);
end

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
% Rename pdf.
movefile("all.pdf", sprintf("%s.pdf", parameters.pdfname));

% Restore the path to oldpath.
setpath(oldpath);  
% Go back to the original directory.
cd(old_dir);

end
