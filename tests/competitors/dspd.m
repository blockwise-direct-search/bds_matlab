function [xopt, fopt, exitflag, output] = dspd(fun, x0, options)
%DSPD (direct search probabilistic descent) solves unconstrained optimization problems without using derivatives.
%
%   It is supported in MATLAB R2017b or later.
%
%   XOPT = BDS(FUN, X0) returns an approximate minimizer XOPT of the function handle FUN, starting the
%   calculations at X0. FUN must accept input X and return a scalar, which is the function value
%   evaluated at X. X0 should be a vector.
%
%   XOPT = BDS(FUN, X0, OPTIONS) performs the computations with the options in OPTIONS. It should be a
%   structure, with the following fields:
%
%   num_blocks                      Number of blocks.
%   MaxFunctionEvaluations          Maximum of function evaluations.
%   MaxFunctionEvaluations_factor   Factor to define the maximum number of 
%                                   function evaluations as a multiplier
%                                   of the dimension of the problem.
%   expand                          Expanding factor of step size.
%   shrink                          Shrinking factor of step size.
%   reduction_factor                Factor of reduction.
%   StepTolerance                   The tolerance for testing whether the step size is small enough.
%   ftarget                         Target of the function value. If the function value is below target,
%                                   then the algorithm terminates.
%   polling_inner                   Polling strategy of each block.
%   direction_set                   direction set of directions.
%   with_cycling_memory             In the opportunistic case (polling_inner == "opportunistic"),
%                                   with_memory decides whether the cycling strategy memorizes
%                                   the history or not.
%   cycling_inner                   Cycling strategy employed in the opportunistic case.
%   accept_simple_decrease          Whether the algorithm accepts simple decrease or not.
%   Algorithm                       Algorithm of BDS. It can be "cbds", "pbds", "rbds", "ds".
%                                   Use Algorithm not algorithm to have the same name as MATLAB.
%   shuffling_period                A positive integer. This is only used for PBDS, which shuffles the blocks
%                                   every shuffling_period iterations.
%   replacement_delay               An integer between 0 and num_blocks-1. This is only used for RBDS. Suppose that
%                                   replacement_delay is r. If block i is selected at iteration k, then it will
%                                   not be selected at iterations k+1, ..., k+r.
%   seed                            Only used by randomized strategy for reproducibility.
%   output_xhist                    Whether the history of points visited is returned or not.
%   output_alpha_hist               Whether the history of step sizes is returned or not.
%   output_block_hist               Whether the history of blocks visited is returned or not.
%
%   [XOPT, FOPT] = BDS(...) also returns the value of the objective function FUN at the
%   solution XOPT.
%
%   [XOPT, FOPT, EXITFLAG] = BDS(...) returns an EXITFLAG that describes the exit
%   condition. The possible values of EXITFLAG are 0, 1, 2, 3.
%
%   0    The StepTolerance of the step size is reached.
%   1    The maximum number of function evaluations is reached.
%   2    The target of the objective function is reached.
%   3    The maximum number of iterations is reached.
%
%   [XOPT, FOPT, EXITFLAG, OUTPUT] = BDS(...) returns a
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


% Set the default value of debug_flag. If options do not contain debug_flag, then
% debug_flag is set to false.
if isfield(options, "debug_flag")
    debug_flag = options.debug_flag;
else
    debug_flag = false;
end

% Check the inputs of the user when debug_flag is true.
if debug_flag
    verify_preconditions(fun, x0, options);
end

% If FUN is a string, then convert it to a function handle.
if ischarstr(fun)
    fun = str2func(fun);
end

% Transpose x0 if it is a row.
x0 = double(x0(:));

% Set the polling directions in D.
n = length(x0);

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
if isfield(options, "num_random_vectors")
    num_directions = max(options.num_random_vectors, 1 + floor(log2(1-log(shrink)/log(expand))));
    choice = true;
else
    choice = false;
end

% Set MAXFUN to the maximum number of function evaluations.
if isfield(options, "MaxFunctionEvaluations")
    MaxFunctionEvaluations = options.MaxFunctionEvaluations;
else
    MaxFunctionEvaluations = get_default_constant("MaxFunctionEvaluations_dim_factor")*n;
end

% Each iteration will at least use one function evaluation. We will perform at most MaxFunctionEvaluations iterations.
% In theory, setting the maximum of function evaluations is not needed. But we do it to avoid infinite
% cycling if there is a bug.
maxit = MaxFunctionEvaluations;

% Set the value of reduction factor.
if isfield(options, "reduction_factor")
    reduction_factor = options.reduction_factor;
else
    reduction_factor = get_default_constant("reduction_factor");
end


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
        alpha_hist = NaN(1, maxit);
    catch
        output_alpha_hist = false;
        warning("The size of alpha_hist exceeds the maximum of memory size limit.")
    end
end

if isfield(options, "alpha_init")
    alpha = options.alpha_init;
else
    alpha = get_default_constant("alpha_init");
end

% Initialize the history of function values.
fhist = NaN(1, MaxFunctionEvaluations);

% Initialize the history of points visited.
if isfield(options, "output_xhist")
    output_xhist = options.output_xhist;
else
    output_xhist = get_default_constant("output_xhist");
end

if output_xhist
    try
        xhist = NaN(n, MaxFunctionEvaluations);
    catch
        output_xhist = false;
        warning("xhist will be not included in the output due to the limit of memory.");
    end
end

% Decide whether to print during the computation.
if isfield(options, "iprint")
    iprint = options.iprint;
else
    iprint = get_default_constant("iprint");
end

% To avoid that the users bring some randomized strings.
if ~isfield(options, "seed")
    random_stream = RandStream("mt19937ar", "Seed", "shuffle");
else
    random_stream = RandStream("mt19937ar", "Seed", options.seed);
end

% Initialize the exitflag where the maximum number of iterations is reached.
exitflag = get_exitflag("MAXIT_REACHED");
xbase = x0;
[fbase, fbase_real] = eval_fun(fun, xbase);
% Set the number of function evaluations.
nf = 1;
if output_xhist
    xhist(:, nf) = xbase;
end
fhist(nf) = fbase_real;
xopt = xbase;
fopt = fbase;

% Check whether FTARGET is reached by FVAL. If it is true, then terminate.
if fopt <= ftarget
    information = "FTARGET_REACHED";
    exitflag = get_exitflag(information);

    % FTARGET has been reached at the very first function evaluation.
    % In this case, no further computation should be entertained, and hence,
    % no iteration should be run.
    maxit = 0;
end

for iter = 1:maxit
    % Record the value of alpha_all of the current iteration in alpha_hist.
    if output_alpha_hist
        alpha_hist(iter) = alpha;
    end

    % Generate the direction set whose directions are uniformly distributed on the unit sphere
    % for each iteration when the Algorithm is "dspd".
    if ~choice
        rv = random_stream.randn(n, 1);
        % Normalize rv.
        rv = rv ./ norm(rv);
        D = [rv, -rv];
    else
        D = random_stream.randn(n, num_directions);
        % Normalize D. vecnorm is introduced in MATLAB 2017a for the first time.
        % Here we are dividing a matrix by a row vector using implicit expansion.
        D = D ./ vecnorm(D);
    end

    suboptions.FunctionEvaluations_exhausted = nf;
    suboptions.MaxFunctionEvaluations = MaxFunctionEvaluations - nf;
    suboptions.cycling_inner = cycling_inner;
    suboptions.with_cycling_memory = with_cycling_memory;
    suboptions.reduction_factor = reduction_factor;
    suboptions.forcing_function = forcing_function;
    suboptions.ftarget = ftarget;
    suboptions.polling_inner = options.polling_inner;
    suboptions.iprint = iprint;

    [sub_xopt, sub_fopt, sub_exitflag, sub_output] = inner_direct_search(fun, xbase,...
        fbase, D, 1:size(D, 2), alpha, suboptions);

    % Update the history of step size.
    if output_alpha_hist
        alpha_hist(:, iter) = alpha_all;
    end

    % Store the history of the evaluations by inner_direct_search,
    % and accumulate the number of function evaluations.
    fhist((nf+1):(nf+sub_output.nf)) = sub_output.fhist;
    if output_xhist
        xhist(:, (nf+1):(nf+sub_output.nf)) = sub_output.xhist;
    end

    nf = nf+sub_output.nf;

    % Update the step sizes and store the history of step sizes.
    % fbase and xbase are used for the computation in the block. fval and
    % xval are always the best function value and point so far.
    if sub_fopt + reduction_factor(3) * forcing_function(alpha) < fbase
        alpha = expand * alpha;
    elseif sub_fopt + reduction_factor(2) * forcing_function(alpha) >= fbase
        alpha = shrink * alpha;
    end

    if (reduction_factor(1) <= 0 && sub_fopt < fbase) || sub_fopt + reduction_factor(1) * forcing_function(alpha) < fbase
        xbase = sub_xopt;
        fbase = sub_fopt;
    end

    % Update xopt and fopt.
    if sub_fopt < fopt
        xopt = sub_xopt;
        fopt = sub_fopt;
    end

    % If sub_output.terminate is true, then inner_direct_search returns
    % boolean value of terminate because either the maximum number of function
    % evaluations or the target of the objective function value is reached.
    % In both cases, the exitflag is set by inner_direct_search.
    if sub_output.terminate
        exitflag = sub_exitflag;
        break;
    end

    % Terminate the computations if the largest component of step size is below a
    % given StepTolerance.
    if alpha < alpha_tol
        exitflag = get_exitflag("SMALL_ALPHA");
        break;
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
    verify_postconditions(fun, xopt, fopt, exitflag, output);
end
