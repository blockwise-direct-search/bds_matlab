function [output] = testbds(parameters)
% output: performance profile;
% parameters: name; solvers_legend; memory; cycling; polling_inner; polling_outer; nb_generator;
% maxfun_dim; maxfun; problems_type; problems_mindim; problems_maxidim; tau;
% clear
% clc

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
% The difference between solvers_legend and name is that solvers_legend must be
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
options_test.is_noisy = parameters.is_noisy;
options_test.noise_level = parameters.noise_level;
% relative: (1+noise_level*noise)*f; absolute: f+noise_level*noise
options_test.noise_abs = parameters.noise_abs;
options_test.noise_type = parameters.noise_type;

options_test.scaling_matrix = false;
options_test.scaling_matrix_factor = 5;

options_solvers.blocks_strategy = parameters.blocks_strategy;
options_solvers.cycling_inner = parameters.cycling_inner;
options_solvers.polling_inner = parameters.polling_inner;
options_solvers.solvers_legend = parameters.solvers_legend;
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
options.path = parameters.path_tests;
path_testdata = fullfile(options.path, 'testdata');
path_outdir = fullfile(options.path, 'testdata', tst);

% mkdir a new folder to save numerical results and source code.
mkdir(path_testdata, tst);
mkdir(path_outdir, 'perf');
options.outdir = fullfile(path_outdir, 'perf');
mkdir(path_outdir, 'src');
path_src = fullfile(path_outdir, 'src');
mkdir(path_outdir, 'tests');
path_tests = fullfile(path_outdir, 'tests');

% Copy the source code and test code to path_outdir.
copyfile(fullfile(parameters.path_src, '*'), path_src);
parameters.path_competitors = fullfile(parameters.path_tests, 'competitors');
copyfile(fullfile(parameters.path_competitors, '*'), path_tests);
copyfile(fullfile(parameters.path_tests, 'private', '*'), path_tests);

% Define current path and target path.
currentPath = parameters.path_tests;
targetPath = path_tests; 

% Get all files under currentPath.
files = dir(fullfile(currentPath, '*'));
 
% Loop all files
for i = 1:length(files)
    if ~files(i).isdir % If it is files, not folders.
        % Construct path of source files and target files
        sourceFile = fullfile(currentPath, files(i).name);
        targetFile = fullfile(targetPath, files(i).name);
        
        % Duplicate files
        copyfile(sourceFile, targetFile);
        
        % Move files
        movefile(targetFile, currentPath);
    end
end

% performance profile
tau = parameters.tau; % Tolerance of convergence test in performance profile
tau_length = length(tau);
options_perf.outdir = options.outdir;
options_perf.stamp = time;
options_perf.solvers = parameters.solvers_legend;
for l = 1:tau_length
    options_perf.tau = tau(l);
    output = perfprof(frec, fmin, options_perf);
end


cd(options.outdir)

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

end

