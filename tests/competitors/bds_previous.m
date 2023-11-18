function [xval, fval, exitflag, output] = bds_previous(fun, x0, options)
%BDS (blockwise direct search) solves unconstrained optimization problems without using derivatives.
%
%   It is supported in MATLAB 2017b or later versions.
%
%   XVAL = BDS(FUN, X0) returns an approximate minimizer XVAL of the function handle FUN, starting the
%   calculations at X0. FUN must accept input X and return a scalar, which is the function value
%   evaluated at X. X0 should be a vector.
%
%   XVAL = BDS(FUN, X0, OPTIONS) performs the computations with the options in OPTIONS. It should be a
%   structure, with the following fields:
%
%   num_blocks                          Number of blocks.
%   maxfun                      Maximum of function evaluations.
%   maxfun_factor               Factor to define the maximum number of function evaluations as a multiplier
%                               of the dimension of the problem.
%   expand                      Expanding factor of step size.
%   shrink                      Shrinking factor of step size.
%   sufficient_decrease_factor  Factor of sufficient decrease condition. Sufficient_decrease_factor(1) is
%                               for the update of fval and xval. Sufficient_decrease_factor(2) and
%                               sufficient_decrease_factor(3) is for the update of step size.
%   StepTolerance               The tolerance for testing whether the step size is small enough.
%   ftarget                     Target of the function value. If the function value is below target,
%                               then the algorithm terminates.
%   polling_inner               Polling strategy of each block.
%   searching_set               Searching set of directions.
%   with_cycling_memory         In the opportunistic case (polling_inner == "opportunistic"),
%                               with_memory decides whether the cycling strategy memorizes
%                               the history or not.
%   cycling_inner               Cycling strategy employed in the opportunistic case.
%   forcing_function            Type of forcing function.
%   Algorithm                   Algorithm of BDS. It can be "cbds", "pbds", "rbds", "ds".
%                               Use Algorithm not algorithm to have the same name as MATLAB.
%   shuffling_period            A positive integer. This is only used for PBDS, which shuffles the blocks
%                               every shuffling_period iterations.
%   replacement_delay           An integer between 0 and num_blocks-1. This is only used for RBDS. Suppose that
%                               replacement_delay is r. If block i is selected at iteration k, then it will
%                               not be selected at iterations k+1, ..., k+r.
%   seed                        Only used by randomized strategy for reproducibility.
%   output_xhist                Whether the history of points visited is returned or not.
%   output_alpha_hist           Whether the history of step sizes is returned or not.
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
%   alpha_hist   History of step size for every iteration.
%   blocks_hist  History of blocks visited.
%   funcCount    The number of function evaluations.
%   message      The information of EXITFLAG.
%
%   Copyright 2023, Haitian Li and Zaikun Zhang.
%   All rights reserved.
%

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

if ~isfield(options, "Algorithm")
    options.Algorithm = get_default_constant("Algorithm");
end

% Get the searching set of directions.
D = get_searching_set(n, options);

% Set the value of expanding factor.
if isfield(options, "expand")
    expand = options.expand;
else
    expand = get_default_constant("expand");
end

% Set the value of shrinking factor.
if isfield(options, "shrink")
    shrink = options.shrink;
else
    shrink = get_default_constant("shrink");
end

% Get the number of directions.
num_directions = size(D, 2);

% Get the number of blocks.
if isfield(options, "num_blocks")
    % The number of directions should be greater or equal to the number of blocks.
    num_blocks = min(num_directions, options.num_blocks);
elseif strcmpi(options.Algorithm, "cbds") || strcmpi(options.Algorithm, "pbds")...
        || strcmpi(options.Algorithm, "rbds")
    % Default value is set as n, which is good for canonical with 2n directions. For
    % other situations, other value may be good.
    num_blocks = n;
elseif strcmpi(options.Algorithm, "ds")
    num_blocks = 1;
end

% Set indices of blocks as 1:num_blocks.
block_indices = 1:num_blocks;

% Set MAXFUN to the maximum number of function evaluations.
if isfield(options, "maxfun_factor") && isfield(options, "maxfun")
    maxfun = min(options.maxfun_factor*n, options.maxfun);
elseif isfield(options, "maxfun_factor")
    maxfun = options.maxfun_factor*n;
elseif isfield(options, "maxfun")
    maxfun = options.maxfun;
else
    maxfun = min(get_default_constant("maxfun"), get_default_constant("maxfun_factor")*n);
end

% Each iteration will at least use one function evaluation. We will perform at most maxfun iterations.
% In theory, setting the maximum of function evaluations is not needed. But we do it to avoid infinite
% cycling if there is a bug.
maxit = maxfun;

% % Set the value of sufficient decrease factor.
% if isfield(options, "sufficient_decrease_factor_level")
%     switch lower(options.sufficient_decrease_factor_level)
%         case "zero"
%             options.sufficient_decrease_factor = [0, 0, 0];
%         case "negligible"
%             options.sufficient_decrease_factor = [1e-16, 1e-16, 1e-16];
%         case "low"
%             options.sufficient_decrease_factor = [1e-8, 1e-8, 1e-8];
%         case "medium"
%             options.sufficient_decrease_factor = [1e-3, 1e-3, 1e-3];
%         case "high"
%             options.sufficient_decrease_factor = [1, 1, 1];
%         case "excessive"
%             options.sufficient_decrease_factor = [10, 10, 10];
%         otherwise
%             error("Unknown sufficient decrease factor level %s", ...
%                 options.sufficient_decrease_factor_level);
%     end
% else
%     if isfield(options, "sufficient_decrease_factor")
%         sufficient_decrease_factor = options.sufficient_decrease_factor;
%     else
%         sufficient_decrease_factor = get_default_constant("sufficient_decrease_factor");
%     end
% end
sufficient_decrease_factor = [0, eps, eps];

% Set the forcing function, which is a function handle.
if isfield(options, "forcing_function")
    forcing_function = options.forcing_function;
else
    forcing_function = get_default_constant("forcing_function");
end

if isfield(options, "forcing_function_type")
    switch options.forcing_function_type
        case "quadratic"
            forcing_function = @(x)x.^2;
        case "cubic"
            forcing_function = @(x)x.^3;
    end
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
    replacement_delay = min(options.replacement_delay, num_blocks-1);
else
    replacement_delay = min(get_default_constant("replacement_delay"), num_blocks-1);
end

% Set the boolean value of WITH_CYCLING_MEMORY.
% WITH_CYCLING_MEMORY is only used when we need to permute the direction_indices. If
% WITH_CYCLING_MEMORY is true, then we will permute the direction_indices by using the
% direction_indices of the previous iteration. Otherwise, we will permute the direction_indices
% with the initial direction_indices of ascending orders.
if isfield(options, "with_cycling_memory")
    with_cycling_memory = options.with_cycling_memory;
else
    with_cycling_memory = get_default_constant("with_cycling_memory");
end

% Initialize the step sizes and alpha_hist, which is the history of step sizes.
if isfield(options, "output_alpha_hist")
    output_alpha_hist = options.output_alpha_hist;
else
    output_alpha_hist = get_default_constant("output_alpha_hist");
end

if output_alpha_hist
    try
        alpha_hist = NaN(num_blocks, maxit);
    catch
        output_alpha_hist = false;
        warning("The size of alpha_hist exceeds the maximum of memory size limit.")
    end
end

if isfield(options, "alpha_init")
    if length(options.alpha_init) == 1
        alpha_all = options.alpha_init*ones(num_blocks, 1);
    elseif length(options.alpha_init) == num_blocks
        alpha_all = options.alpha_init;
    else
        error("The length of alpha_init should be equal to num_blocks or equal to 1.");
    end
    % Try alpha_all = 0.5 * max(abs(x0), 1) in the canonical case.
elseif isfield(options, "alpha_init_perturbed") && options.alpha_init_perturbed
    alpha_all = 0.00025 * ones(num_blocks, 1);
    alpha_all(x0 ~= 0) = 1.05 * abs(x0(x0 ~= 0));
    % alpha_all = 0.5 * max(abs(x0), ones(num_blocks, 1));
else
    alpha_all = ones(num_blocks, 1);
end

% Decide which polling direction belongs to which block.
searching_set_indices = divide_searching_set(num_directions, num_blocks);

% Initialize the history of function values.
fhist = NaN(1, maxfun);

% Initialize the history of points visited.
if isfield(options, "output_xhist")
    output_xhist = options.output_xhist;
else
    output_xhist = get_default_constant("output_xhist");
end

if output_xhist
    try
        xhist = NaN(n, maxfun);
    catch
        output_xhist = false;
        warning("xhist will be not included in the output due to the limit of memory.");
    end
end

if isfield(options, "output_block_hist")
    output_block_hist = options.output_block_hist;
else
    output_block_hist = get_default_constant("output_block_hist");
end

% Initialize the history of blocks visited.
block_hist = NaN(1, maxfun);
xval = x0;
fval = eval_fun(fun, xval);
% Set the number of function evaluations.
nf = 1;
if output_xhist
    xhist(:, nf) = xval;
end
fhist(nf) = fval;

% Check whether FTARGET is reached by FVAL. If it is true, then terminate.
if fval <= ftarget
    information = "FTARGET_REACHED";
    exitflag = get_exitflag(information);

    % FTARGET has been reached at the very first function evaluation.
    % In this case, no further computation should be entertained, and hence,
    % no iteration should be run.
    maxit = 0;
end

% To avoid that the users bring some randomized strings.
if ~isfield(options, "seed")
    random_stream = RandStream("mt19937ar", "Seed", "shuffle");
else
    random_stream = RandStream("mt19937ar", "Seed", options.seed);
end


% Start the actual computations.
for iter = 1:maxit
    % Record the value of alpha_all of the current iteration in alpha_hist.
    if output_alpha_hist
        alpha_hist(:, iter) = alpha_all;
    end

    % Shuffle the blocks every shuffling_period iterations.
    % Why iter-1? Since we will permute block_indices at the initial stage.
    if strcmpi(options.Algorithm, "pbds") && mod(iter - 1, shuffling_period) == 0
        % Make sure that shuffling_period is defined when the Algorithm is "pbds".
        block_indices = random_stream.randperm(num_blocks);
    end

    % Get the block that is going to be visited.
    if strcmpi(options.Algorithm, "rbds")
        % If replacement_delay is 0, then select a block randomly from block_indices for
        % each iteration. If iter is equal to 1, then the block that we are going to visit
        % is selected randomly from block_indices.
        if replacement_delay == 0 || iter == 1
            block_indices = random_stream.randi([1, num_blocks]);
        else
            % Record the number of blocks visited.
            num_visited = sum(~isnan(block_hist));
            % Get the number of blocks that we are going to exclude in the following selection.
            block_visited_slices_length = min(num_visited, replacement_delay);
            % Get the indices of blocks that we are going to exclude in the following selection.
            block_visited_slices = block_hist(num_visited-block_visited_slices_length+1:num_visited);
            % Set the default value of initial block_indices.
            block_initial_indices = 1:num_blocks;
            % Remove elements of block_indices appearing in block_visited_slice.
            block_real_indices = block_initial_indices(~ismember(block_initial_indices, block_visited_slices));
            % Generate a random index from block_real_indices.
            idx = random_stream.randi(length(block_real_indices));
            block_indices = block_real_indices(idx);
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
        suboptions.forcing_function = forcing_function;
        suboptions.ftarget = ftarget;
        suboptions.polling_inner = options.polling_inner;

        [xval, fval, sub_exitflag, suboutput] = inner_direct_search_copy(fun, xval,...
            fval, D(:, direction_indices), direction_indices,...
            alpha_all(i_real), suboptions);

        % Update the history of step size.
        if output_alpha_hist
            alpha_hist(:, iter) = alpha_all;
        end

        % Store the history of the evaluations by inner_direct_search,
        % and accumulate the number of function evaluations.
        fhist((nf+1):(nf+suboutput.nf)) = suboutput.fhist;
        if output_xhist
            xhist(:, (nf+1):(nf+suboutput.nf)) = suboutput.xhist;
        end
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

        % Update the step sizes and store the history of step sizes.
        success = suboutput.success;
        reduction = suboutput.reduction;
        if success
            alpha_all(i_real) = expand * alpha_all(i_real);
        else
            if ~reduction
                alpha_all(i_real) = shrink * alpha_all(i_real);
            end
        end

        % Terminate the computations if the largest component of step size is below a
        % given StepTolerance.
        if max(alpha_all) < alpha_tol
            terminate = true;
            exitflag = get_exitflag("SMALL_ALPHA");
            break
        end

        % Retrieve the order of the polling directions in the searching set.
        searching_set_indices{i_real} = suboutput.direction_indices;

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

% Truncate HISTORY into a vector of nf length.
output.funcCount = nf;
output.fhist = fhist(1:nf);
if output_xhist
    output.xhist = xhist(:, 1:nf);
end
if output_alpha_hist
    output.alpha_hist = alpha_hist(1:min(iter, maxit));
end

% Record the number of blocks visited.
num_blocks_visited = sum(~isnan(block_hist));

% Record the blocks visited.
if output_block_hist
    output.blocks_hist = block_hist(1:num_blocks_visited);
end


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

% verify_postconditions is to detect whether the output is in the right form when debug_flag is true.
if debug_flag
    verify_postconditions(fun, xval, fval, exitflag, output);
end