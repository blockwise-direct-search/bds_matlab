function set_StepTolerance(Algorithm, options)
% Get the ratio of StepTolerance and the minimum of gradient value among
% fhist.
%

if nargin < 2
    options = struct();
end
options.solver_name = "bds";
options.Algorithm = Algorithm;

% Record the current path.
oldpath = path();
% Restore the "right out of the box" path of MATLAB.
restoredefaultpath;
% Record the current directory.
old_dir = pwd();

exception = [];

try

    % Add the paths that we need to use in the performance profile into the MATLAB
    % search path.
    current_path = mfilename("fullpath");
    path_tests = fileparts(current_path);
    path_root = fileparts(path_tests);
    path_src = fullfile(path_root, "src");
    addpath(path_root);
    addpath(path_tests);
    addpath(path_src);

    % If the folder of testdata does not exist, make a new one.
    path_testdata = fullfile(path_tests, "testdata");
    if ~exist(path_testdata, "dir")
        mkdir(path_testdata);
    end

    % Tell MATLAB where to find MatCUTEst.
    locate_matcutest();

    % Get list of problems
    s.type = 'u'; % Unconstrained: 'u'
    s.mindim = 1; % Minimum of dimension
    s.maxdim = 100; % Maximum of dimension
    s.blacklist = [];
    s.blacklist = [s.blacklist, { 'ARGTRIGLS', 'BROWNAL', ...
        'COATING', 'DIAMON2DLS', 'DIAMON3DLS', 'DMN15102LS', ...
        'DMN15103LS', 'DMN15332LS', 'DMN15333LS', 'DMN37142LS', ...
        'DMN37143LS', 'ERRINRSM', 'HYDC20LS', 'LRA9A', ...
        'LRCOVTYPE', 'LUKSAN12LS', 'LUKSAN14LS', 'LUKSAN17LS', 'LUKSAN21LS', ...
        'LUKSAN22LS', 'MANCINO', 'PENALTY2', 'PENALTY3', 'VARDIM',
        }];
   

    problem_names = secup(s);
    num_problems = length(problem_names);
    data = cell(num_problems, 4);

    % Set output_xhist to be true to have output.xhist for calculating the
    % norm of the gradient.
    options.output_xhist = true;

    % Set maxfun large enough to see the ratio between the norm of the
    % gradient and the StepTolerance.
    if ~isfield(options, "maxfun")
        options.maxfun = 1e5;
    end
    
    % Use time to distinguish.
    time_str = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm'));
    % Trim time string.
    time_str = trim_time(time_str);
    % Rename tst as mixture of time Algorithm and StepTolerance.
    tst = strcat("ratio", "_", time_str);
    path_testdata = fullfile(path_tests, "testdata");
    path_ratio = fullfile(path_testdata, tst);

    % Make a new folder to save the results and source code.
    mkdir(path_ratio);
    path_ratio_src = fullfile(path_ratio, "src");
    mkdir(path_ratio_src);
    copyfile(fullfile(path_src, "*"), path_ratio_src);

    % Make a txt file to store the ratio that are recorded.
    problem_names_str = strings(1, length(problem_names));
    for i = 1:num_problems
        problem_names_str(i) = problem_names{i};
    end
    filePath = strcat(path_ratio, "/ratio.txt");
    fileID = fopen(filePath, 'w');
    fprintf(fileID, '%s\n', options.Algorithm);
    fprintf(fileID, '%s\n', num2str(options.maxfun));
    if isfield(options, "StepTolerance")
        fprintf(fileID, '%s\n', num2str(options.StepTolerance));
    end
    for i_problem = 1:num_problems
        p = macup(problem_names(1, i_problem));
        ratio = test_gradient_stepsize(p.name, options);
        fprintf(fileID, '%s %s\n', p.name, num2str(ratio));
    end
    fclose(fileID);

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