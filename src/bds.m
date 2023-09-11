function [xval, fval, exitflag, output] = bds(fun, x0, options)
%BDS (blockwise direct search) solves unconstrained optimization problems without using derivatives. 
%
%   TODO: here(MATLAB version, what is bds)
%
%   XVAL = BDS(FUN, X0) returns an approximate minimizer XVAL of the function handle FUN, starting the
%   calculations at X0. FUN must accept input X and returns a scalar, which is the function value
%   evaluated at X. X0 should be a vector.
%
%   XVAL = BDS(FUN, X0, OPTIONS) performs the computations with the options in OPTIONS. It should be a
%   structure, with the following fields:
%   
%   nb                          Number of blocks.
%   maxfun                      Maximum of function evaluations.
%   maxfun_dim                  Factor to define maximum number of function evaluations as a multiplier
%                               of the dimension of the problem.    
%   expand                      Expanding factor of step size.
%   shrink                      Shrinking factor of step size.
%   sufficient_decrease_factor  Factor of sufficient decrease condition.
%   StepTolerance               The tolerance for testing whether the step size is small enough.
%   ftarget                     Target of function value. If function value is below ftarget, 
%                               then the algorithm terminates.
%   polling_inner               Polling strategy of each block.
%   with_cycling_memory         In the opportunistic case (polling_inner == "opportunistic"), 
%                               with_meory decides whether the cycling strategy memorizes 
%                               the history or not.
%   cycling_inner               Cycling strategy employed in the opportunistic case.
%   accept_simple_decrease      Whether the algorithm accepts simple decrease or not.
%   Algorithm                   Algorithm of BDS. It can be "cbds", "pbds", "rbds", "dspd", "ds".
%                               Use Algorithm not algorithm to have the same name as MATLAB.
%   shuffling_period            A positive integer. This is only used for PBDS, which shuffles the blocks
%                               every shuffling_period iterations.    
%   replacement_delay           An integer between 0 and nb-1. This is only used for RBDS. Suppose that 
%                               replacement_delay is r. If block i is selected at iteration k, then it will 
%                               not be selected at iterations k+1, ..., k+r. 
%   seed                        Only used by randomized strategy for reproducibility.
%
%   [XVAL, FVAL] = BDS(...) also returns the value of the objective function FUN at the 
%   solution XVAL.
%
%   [XVAL, FVAL, EXITFLAG] = BDS(...) returns an EXITFLAG that describes the exit 
%   condition. The possible values of EXITFLAG are 0, 1, 2, 3.
%
%   0    The StepTolerance of the step size is reached.
%   1    The maximum number of function evaluations is reached.
%   2    The target of the objective function is reached.
%   3    The maximum number of iterations is reached.
%
%   [XVAL, FVAL, EXITFLAG, OUTPUT] = BDS(...) returns a
%   structure OUTPUT with the following fields: 
%
%   fhist        History of function values.
%   xhist        History of points visited.
%   alpha_hist   History of step size.
%   blocks_hist  History of blocks visited.
%   funcCount    The number of function evaluations.
%   message      The information of EXITFLAG.
%
%   Copyright 2023, Haitian Li and Zaikun Zhang.
%   All rights reserved.
%

% Get the information of current MATLAB.
matlabVersion = ver('MATLAB');

% Extract the main version and the minor version of MATLAB.
versionParts = strsplit(matlabVersion.Version, '.');

% Change the version of matlab into floating number.
majorVersion = str2double(versionParts{1});
minorVersion = str2double(versionParts{2});

% Check whether the verion of MATLAB is released lower than 2017a.
if majorVersion < 9 || (majorVersion == 9 && minorVersion < 1)
    error("The verison of MATLAB is low, please update to 2017a or later versions")
end

% Set options to an empty structure if it is not provided.
if nargin < 3
    options = struct();
end

% Check the inputs of the user when debug_flag is true.
debug_flag = is_debugging();
if debug_flag
    verify_preconditions(fun, x0, options);
end

% If FUN is a string, then convert it to a function handle.
if ischarstr(fun)
    fun = str2func(fun);
end

% We set the initial flag to NaN. This value will be modified by procedures.
% If EXITFLAG is set to NaN on exit, it means that there is a bug.
exitflag = NaN;

% Transpose x0 if it is a row.
x0 = double(x0(:));

% Set the polling directions in D.
n = length(x0);

if isfield(options, "seed")
    seed = options.seed;
else
    seed = get_default_constant("seed");
end

if ~isfield(options, "Algorithm")
    options.Algorithm = get_default_constant("Algorithm");
end

if strcmpi(options.Algorithm, "cbds") || strcmpi(options.Algorithm, "pbds")...
        || strcmpi(options.Algorithm, "ds") || strcmpi(options.Algorithm, "rbds")
    D = get_searching_set(n, options);
end

% Set the value of expand factor.
if isfield(options, "expand")
    expand = options.expand;
else
    expand = get_default_constant("expand");
end

% Set the value of shrink factor.
if isfield(options, "shrink")
    shrink = options.shrink;
else
    shrink = get_default_constant("shrink");
end

% Get the number of directions.
if strcmpi(options.Algorithm, "dspd")
    if isfield(options, "num_random_vectors")
        m = max(options.num_random_vectors, ceil(log2(1-log(shrink))/log(expand)));
    else
        m = max(get_default_constant("num_random_vectors"), ceil(log2(1-log(shrink))/log(expand)));
    end
else
    m = size(D, 2);
end
 
% Get the number of blocks.
if isfield(options, "nb")
    nb = options.nb;
elseif strcmpi(options.Algorithm, "cbds") || strcmpi(options.Algorithm, "pbds")...
        || strcmpi(options.Algorithm, "rbds")
    % Default value is set as n, which is good for canonical with 2n directions. For
    % other situations, other value may be good.
    nb = n;
elseif strcmpi(options.Algorithm, "dspd") || strcmpi(options.Algorithm, "ds")
    nb = 1;
end

% The number of directions should be greater or equal to the number of blocks.
nb = min(m, nb);

% Set indices of blocks as 1:nb.
block_indices = 1:nb;

% Set MAXFUN to the maximum number of function evaluations.
if isfield(options, "maxfun_dim") && isfield(options, "maxfun")
    maxfun = min(options.maxfun_dim*n, options.maxfun);
elseif isfield(options, "maxfun_dim")
    maxfun = options.maxfun_dim*n;
elseif isfield(options, "maxfun")
    maxfun = options.maxfun;
else
    maxfun = min(get_default_constant("maxfun"), get_default_constant("maxfun_dim")*n);
end

% Each iteration will at least use one function evaluation. We will perform at most maxfun iterations.
% In theory, setting the maximum of function evaluations is not needed. But we do it to avoid infinite 
% cyling if there is a bug.
maxit = maxfun;

% Set the value of sufficient decrease factor.
if isfield(options, "sufficient_decrease_factor")
    sufficient_decrease_factor = options.sufficient_decrease_factor;
else
    sufficient_decrease_factor = get_default_constant("sufficient_decrease_factor");
end

% Set the boolean value of accept_simple_decrease, which is for updating xval and fval, but not 
% for stepsize. 
if isfield(options, "accept_simple_decrease")
    accept_simple_decrease = options.accept_simple_decrease;
else
    accept_simple_decrease = get_default_constant("accept_simple_decrease");
end

% Set the value of StepTolerance. The algorithm will terminate if the stepsize is less than 
% the StepTolerance.
if isfield(options, "StepTolerance")
    alpha_tol = options.StepTolerance;
else
    alpha_tol = get_default_constant("StepTolerance");
end

% Set the target of the objective function.
if isfield(options, "ftarget")
    ftarget = options.ftarget;
else
    ftarget = get_default_constant("ftarget");
end

% Set the value of polling_inner. This is the polling strategy employed within one block.
if ~isfield(options, "polling_inner")
    options.polling_inner = get_default_constant("polling_inner");
end

% Set the value of cycling_inner, which represents the cycling strategy inside each block.
if isfield(options, "cycling_inner")
    cycling_inner = options.cycling_inner;
else
    cycling_inner = get_default_constant("cycling_inner");
end

% Set the value of shuffling_period, which shuffles the blocks every shuffling_period iterations.
if strcmpi(options.Algorithm, "pbds") && isfield(options, "shuffling_period")
    shuffling_period = options.shuffling_period;
else
    shuffling_period = get_default_constant("shuffle_period");
end

% Set the value of replacement_delay. The default value of replacement_delay is set to 0. 
if strcmpi(options.Algorithm, "rbds") && isfield(options, "replacement_delay")
    replacement_delay = min(options.replacement_delay, nb-1);
else
    replacement_delay = min(get_default_constant("replacement_delay"), nb-1);
end

% Set the boolean value of WITH_CYCLING_MEMORY. 
% TODO: write it more detailedly.
if isfield(options, "with_cycling_memory")
    with_cycling_memory = options.with_cycling_memory;
else
    with_cycling_memory = get_default_constant("with_cycling_memory");
end

% Initialize the step sizes and alpha_hist, which is the history of step sizes.
alpha_hist = NaN(nb, maxit);
if isfield(options, "alpha_init")
    alpha_all = options.alpha_init*ones(nb, 1);
else
    alpha_all = ones(nb, 1);
end

% Decide which polling direction belongs to which block.
searching_set_indices = divide_searching_set(m, nb);

% Initialize the history of function values.
fhist = NaN(1, maxfun);

% Initialize the history of points visited.
xhist = NaN(n, maxfun); 

% Initialize the history of blocks visited.
block_hist = NaN(1, maxfun);
xval = x0; 
fval = eval_fun(fun, xval);
% Set the number of function evaluations.
nf = 1; 
fhist(nf) = fval;
xhist(:, nf) = xval;

% Check whether FTARGET is reached by FVAL. If it is true, then terminate.
if fval <= ftarget
    information = "FTARGET_REACHED";
    exitflag = get_exitflag(information);
    
    % FTARGET has been reached at the very first function evaluation. 
    % In this case, no further computation should be entertained, and hence, 
    % no iteration should be run.
    maxit = 0;
end

% To avoid that the users bring some randomized string.
% TODO: mention objective function and random generator
random_stream = rng(seed);

% Start the actual computations.
for iter = 1:maxit
    % Record the value of alpha_all of the current iteration in alpha_hist.
    alpha_hist(:, iter) = alpha_all;
    
    % Shuffle the blocks every shuffling_period iterations.
    % Why iter-1? Since we will permute block_indices at the initial stage.
    if strcmpi(options.Algorithm, "pbds") && mod(iter - 1, shuffling_period) == 0
        % Make sure that shuffling_period is defined when Algorithm is "sbds".
        block_indices = randperm(nb);
    end
    
    % Get the block that is going to be visited.
    if strcmpi(options.Algorithm, "rbds")
        % If replacement_delay is 0, then select a block randomly from block_indices for 
        % each iteration. If iter is equal to 1, then the block that we are going to visit
        % is selected randomly from block_indices.
        if replacement_delay == 0 || iter == 1
            block_indices = randi([1, nb]);
        else
            % Record the number of blocks visited.
            num_visited = sum(~isnan(block_hist));
            % Get the number of blocks that we are going to exclude in the following selection.
            block_visited_slices_length = min(num_visited, replacement_delay);
            % Get the indices of blocks that we are going to exclude in the following selection.
            block_visited_slices = block_hist(num_visited-block_visited_slices_length+1:num_visited);
            % Set default value of initial block_indices.
            block_initial_indices = 1:nb;
            % Remove elements of block_indices appearing in block_visited_slice.
            block_real_indices = block_initial_indices(~ismember(block_initial_indices, block_visited_slices));
            % Generate a random index from block_real_indices.
            idx = randi(length(block_real_indices));
            block_indices = block_real_indices(idx);
        end
    end
    
    % Generate the searching set whose directions are uniformly distributed on the unit sphere
    % for each iteration when options.Algorithm is "dspd".
    if strcmpi(options.Algorithm, "dspd")
        if m == 2
            rv = randn(n, 1);
            % Normalize rv.
            rv = rv ./ norm(rv);
            D = [rv, -rv];
        else
            D = randn(m, n);
            % Normalize D. vecnorm is introduced in MATLAB 2017a for the first time.
            % Here we are dividing a matrix by a row vector using implicit expansion.
            D = D ./ vecnorm(D);
        end
    end

    for i = 1:length(block_indices)
        % If block_indices is 1 3 2, then block_indices(2) = 3, which is the real block that we are
        % going to visit.
        i_real = block_indices(i);
        
        % Record the number of blocks visited.
        num_visited = sum(~isnan(block_hist));

        % Record the block that is going to be visited.
        block_hist(num_visited+1) = i_real;
        
        % Get indices of directions in the i-th block.
        direction_indices = searching_set_indices{i_real}; 
        
        suboptions.maxfun = maxfun - nf;
        suboptions.cycling = cycling_inner;
        suboptions.with_cycling_memory = with_cycling_memory;
        suboptions.sufficient_decrease_factor = sufficient_decrease_factor;
        suboptions.ftarget = ftarget;
        suboptions.polling_inner = options.polling_inner;
        suboptions.accept_simple_decrease = accept_simple_decrease;
        
        [xval, fval, sub_exitflag, suboutput] = inner_direct_search(fun, xval,...
            fval, D(:, direction_indices), direction_indices,...
            alpha_all(i_real), suboptions);
        
        % Update the history of step size.
        alpha_hist(:, iter) = alpha_all;
        
        % Store the history of the evaluations by inner_direct_search, 
        % and accumulate the number of function evaluations.
        fhist((nf+1):(nf+suboutput.nf)) = suboutput.fhist;
        xhist(:, (nf+1):(nf+suboutput.nf)) = suboutput.xhist;
        nf = nf+suboutput.nf;
        
        % If suboutput.terminate is true, then inner_direct_search returns 
        % boolean value of terminate because either the maximum number of function
        % evaluations or the target of the objective function value is reached. 
        % In both cases, the exitflag is set by inner_direct_search.
        terminate = suboutput.terminate;
        if terminate
            exitflag = sub_exitflag;
            break;
        end
        
        % Retrieve the order of the polling directions and check whether a
        % sufficient decrease has been achieved in inner_direct_search.
        searching_set_indices{i_real} = suboutput.direction_indices;
        success = suboutput.success;
        
        % Update the step sizes and store the history of step sizes.
        if success
            alpha_all(i_real) = expand * alpha_all(i_real);
        else
            alpha_all(i_real) = shrink * alpha_all(i_real);
        end
        
        % Terminate the computations if the largest component of step size is below a
        % given StepTolerance.
        if max(alpha_all) < alpha_tol
            terminate = true;
            exitflag = get_exitflag("SMALL_ALPHA");
            break
        end
    end
    
    % Check whether one of SMALL_ALPHA, MAXFUN_REACHED, and FTARGET_REACHED is reached.
    if terminate
        break;
    end
    
    % Check whether MAXIT is reached.
    if iter == maxit
        exitflag = get_exitflag("MAXIT_REACHED");
    end
    
end

rng(seed);

% Truncate HISTORY into an nf length vector.
output.funcCount = nf;
output.fhist = fhist(1:nf);
output.xhist = xhist(:, 1:nf);
output.alpha_hist = alpha_hist(1:min(iter, maxit));

% Record the number of blocks visited.
num_blocks_visited = sum(~isnan(block_hist));

% Record the blocks visited.
output.blocks_hist = block_hist(1:num_blocks_visited);

switch exitflag
    case {get_exitflag("SMALL_ALPHA")}
        output.message = "The StepTolerance of the step size is reached.";
    case {get_exitflag("MAXFUN_REACHED")}
        output.message = "The maximum number of function evaluations is reached.";
    case {get_exitflag("FTARGET_REACHED")}
        output.message = "The target of the objective function is reached.";
    case {get_exitflag("MAXIT_REACHED")}
        output.message = "The maximum number of iterations is reached.";
    otherwise
        output.message = "Unknown exitflag";
end

% verify_postconditions is to detect whether output is in right form when debug_flag is true.
if debug_flag
    verify_postconditions(fun, xval, fval, exitflag, output);
end
