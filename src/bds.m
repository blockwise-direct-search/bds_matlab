function [xval, fval, exitflag, output] = bds(fun, x0, options)
%   BDS Unconstrained nonlinear minimization (direct search with blocks).
%
%   XVAL = BDS(FUN, X0) starts at X0 and attempts to find
%   local minimizer X of the function FUN.  FUN is a function handle. FUN
%   accepts input X and returns a scalar function value F evaluated at X.
%   X0 should be a vector.
%
%   XVAL = BDS(FUN, X0, OPTIONS) minimizes with the
%   default optimization parameters replaced by values in the structure OPTIONS,
%   BLOCKWISE_DIRECT_SEARCH uses these options: nb, maxfun, maxfun_dim,
%   expand, shrink, sufficient decrease factor, StepTolerance, ftarget, polling_inner,
%   blocks_strategy, with_memory, cycling, accept_simple_decrease.
%
%   nb - number of blocks
%   maxfun - maximum of function evaluation
%   maxfun_dim - factor of maximum of function evaluation regarding to
%               dimenstions.
%   expand - expanding factor of step size
%   shrink - shrinking factor of step size
%   sufficient_decrease_factor - factor of sufficient decrease condition
%   StepTolerance - StepTolerance of step size. If step size is below StepTolerance, then the
%        algorithm terminates.
%   ftarget - If function value is below ftarget, then the algorithm terminates.
%   polling_inner - polling strategy of indices in one block
%
%   [XVAL, FVAL] = BDS(...) returns the value of the
%   objective function, described in FUN, at XVAL.
%
%   [XVAL, FVAL, EXITFLAG] = BDS(...) returns an EXITFLAG
%   that describes the exit condition. The information of EXITFLAG will be
%   given in output.message.
%
%   [XVAL, FVAL, EXITFLAG, OUTPUT] = BDS(...) returns a
%   structure OUTPUT with fields
%
%   fhist      History of function value
%   xhist      History of points that being calculated
%   alpha_hist History of step size every iteration
%   funcCount  The number of function evaluations.
%
%   The number of function evaluations is OUTPUT.funcCount.
%   The history of function evaluation is OUTPUT.fhist.
%   The history of points is OUTPUT.xhist and
%   The l2-norm of gradient in OUTPUT.ghist (only for CUTEst problems).

% Set options to an empty structure if it is not supplied.
if nargin < 3
    options = struct();
end

% Verify_preconditions: If debug_flag is true, then verify_preconditions is to verify
% input. If input_correctness is false, then assert may let the code crash.
debug_flag = is_debugging();
if debug_flag
    verify_preconditions(fun, x0, options);
end

% The exit flag will be set at each possible exit of the algorithm.
% Therefore, if it is set to NaN on exit, it means that there is a bug.
exitflag = NaN;

% Transpose x0 if it is a row.
x0 = double(x0(:));

% Set the polling directions in D.
n = length(x0);
D = get_searching_set(n, options);
m = size(D, 2); % number of directions

% Set the default number of blocks.
if isfield(options, "nb")
    nb = options.nb;
else
    % Default value is good for canonical with 2n directions. For
    % other situations, other value may be good.
    nb = n;
end

% If number of directions is less than number of blocks, then the number of
% blocks is defined as the number of directions.
nb = min(m, nb);
block_indices = 1:nb;

% Set maxfun to the maximum number of function evaluations. The default
% value is 1e5.

if isfield(options, "maxfun_dim") && isfield(options, "maxfun")
    maxfun = min(options.maxfun_dim*n, options.maxfun);
elseif isfield(options, "maxfun_dim")
    maxfun = options.maxfun_dim*n;
elseif isfield(options, "maxfun")
    maxfun = options.maxfun;
else
    maxfun = min(get_default_constant("maxfun"), get_default_constant("maxfun_dim")*n);
end

% Set the maximum of iterations. If complete polling is used, then the
% maximum number of iterations given below CANNOT be reached. If the
% opportunistic case is used, then the maximum number of iterations may
% be reached (although, we hope that it does not).
% ceil(10*maxfun/m) may not be enough. Since there are some cases that
% maxit is exhausted and other terminations are not reached.
maxit = maxfun;

% Set the default expanding factor.
if isfield(options, "expand")
    expand = options.expand;
else
    expand = get_default_constant("expand");
end

% Set the default shrinking factor.
if isfield(options, "shrink")
    shrink = options.shrink;
else
    shrink = get_default_constant("shrink");
end

% Set the default sufficient decrease factor.
if isfield(options, "sufficient_decrease_factor")
    sufficient_decrease_factor = options.sufficient_decrease_factor;
else
    sufficient_decrease_factor = get_default_constant("sufficient_decrease_factor");
end

% Set the default boolean value of accept_simple_decrease. If
% accept_simple_decrease is set to be true, it means the algorithm accepts
% simple decrease to update xval and fval. However, alpha is always updated
% by whether meeting sufficient decrease.
if isfield(options, "accept_simple_decrease")
    accept_simple_decrease = options.accept_simple_decrease;
else
    accept_simple_decrease = get_default_constant("accept_simple_decrease");
end

% Set the default StepTolerance of step size. If the step size reaches a value
% below this StepTolerance, then the algorithm is stopped.
if isfield(options, "StepTolerance")
    alpha_tol = options.StepTolerance;
else
    alpha_tol = get_default_constant("StepTolerance");
end

% Set the target on the objective function. If an evaluation of the
% objective function is below the target (the problem is unconstrained),
% then the algorithm is stopped.
if isfield(options, "ftarget")
   ftarget = options.ftarget;
else
   ftarget = get_default_constant("ftarget");
end

% Set the default blocks_strategy. Default one is Gauss-Seidel.
if ~isfield(options, "blocks_strategy")
    options.blocks_strategy = get_default_constant("blocks_strategy");
end

% Set the default inner polling strategy. This is the polling strategy
% employed within one block.
if ~isfield(options, "polling_inner")
    options.polling_inner = get_default_constant("polling_inner");
end

if isfield(options, "cycling_inner")
    cycling_inner = options.cycling_inner;
else
    cycling_inner = get_default_constant("cycling_inner");
end


% Set the default value for the boolean indicating whether the cycling
% strategy employed in the opportunistic case memorizes the history or not.
if isfield(options, "with_memory")
    with_memory = options.with_memory;
else
    with_memory = get_default_constant("with_memory");
end

% Set initial step size and alpha_hist to store the history of step size.
alpha_hist = NaN(nb, maxit);
if isfield(options, "alpha_init")
    alpha_all = options.alpha_init*ones(nb, 1);
else
    alpha_all = ones(nb, 1);
end

% Divide the indices of the polling directions for each block.
searching_set_indices = divide_searching_set(m, nb);

% Initialize the computations.
fhist = NaN(1, maxfun); % history of function values
xhist = NaN(n, maxfun); % history of points having been visited
hist.block = zeros(1, maxfun); % history of blocks having been visited
xval = x0; % current iterate
fval = fun(xval);
nf = 1; % number of function evaluations
fhist(nf) = fval;
xhist(:, nf) = xval;


% Check whether ftarget is reached by fval. If this is the case, the
% computations afterwards should NOT be done.
if fval <= ftarget
    information = "FTARGET_REACHED";
    exitflag = get_exitflag(information);

    % The target function value has been reached at the very first function
    % evaluation. In this case, no further computation should be
    % entertained, and hence, no iteration should be run.
    maxit = 0;
end

% The number of blocks having been visited. When we store alpha_hist, this
% parameter is needed.
nb_visited = 0;

% Start the actual computations.
% nb blocks have been explored after the number of iteration goes from k to k+1.

for iter = 1 : maxit
    % record the value of alpha_all of the current iteration in alpha_hist.
    alpha_hist(:, iter) = alpha_all;

    % Let xbase be the point from which the polling directions are
    % employed. In one iteration, all the block use the same base point.
    % The corresponding value of the objective function is stored in fbase.
    xbase = xval(:);
    fbase = fval;

    block_indices = permutation(block_indices, options);
    options.permutation_indicator = false;

    for i = 1:nb
        % In case of permutation.
        i_real = block_indices(i);

        direction_indices = searching_set_indices{i_real}; % get indices in the i-th block

        suboptions.maxfun = maxfun - nf;
        % Memory and cycling are needed since we permutate indices in inner_direct_search
        suboptions.cycling = cycling_inner;
        suboptions.with_memory = with_memory;
        suboptions.sufficient_decrease_factor = sufficient_decrease_factor;
        suboptions.ftarget = ftarget;
        suboptions.polling_inner = options.polling_inner;
        suboptions.accept_simple_decrease = accept_simple_decrease;
        
        [xval, fval, sub_exitflag, suboutput] = inner_direct_search(fun, xval,...
            fval, D(:, direction_indices), direction_indices,...
            alpha_all(i_real), suboptions);

        % The i-th block has been visited recently.
        hist.block(iter) = i_real;
        % Update the history of step size.
        alpha_hist(:, iter) = alpha_all;
        % Update the number of blocks having been visited.
        nb_visited = nb_visited + 1;

        % Store the history of the evaluations performed by
        % inner_direct_search, and adjust the number of function
        % evaluations.
        fhist((nf+1):(nf+suboutput.nf)) = suboutput.fhist;
        xhist(:, (nf+1):(nf+suboutput.nf)) = suboutput.xhist;
        nf = nf+suboutput.nf;

        % If suboutput.terminate is true, then inner_direct_search returned
        % because either the maximum number of function evaluations or the
        % target on the objective function value is reached. In both cases,
        % the exitflag is set by inner_direct_search.
        terminate = suboutput.terminate;
        if terminate
            exitflag = sub_exitflag;
            break;
        end

        % Retrieve the order the polling direction and check whether a
        % sufficient decrease has been achieved in inner_direct_search.
        searching_set_indices{i_real} = suboutput.direction_indices;
        success = suboutput.success;

        % Update the step sizes and store the history of step sizes.
        if success
            alpha_all(i_real) = expand * alpha_all(i_real);
        else
            alpha_all(i_real) = shrink * alpha_all(i_real);
        end

        % Terminate the computations if the largest step size is below a
        % given StepTolerance.
        if max(alpha_all) < alpha_tol
            terminate = true;
            exitflag = get_exitflag("SMALL_ALPHA");
            break
        end
    end

    % The following case can be reached (SMALL_ALPHA, MAXFUN_REACHED,
    % FTARGET_REACHED).
    if terminate
        break;
    end

    % Set the exit flag corresponding to "MAXIT_REACHED" on the last
    % iteration. Note that it should be set at last, because another
    % stopping criterion may be reached at the last iteration.
    if iter == maxit
        exitflag = get_exitflag("MAXIT_REACHED");
    end

end

% Set useful pieces on information about the solver's history in output.
output.funcCount = nf;
output.fhist = fhist(1:nf);
output.xhist = xhist(:, 1:nf);
output.alpha_hist = alpha_hist(:, 1:min(iter, maxit));

% Verify_postconditions: If debug_flag is true, then verify_postconditions is
% to verify output. If output_correctness is false, then assert will let code crash.
if debug_flag
    verify_postconditions(fun, xval, fval, exitflag, output);
end

switch exitflag
    case {get_exitflag("SMALL_ALPHA")}
        output.message = "The StepTolerance on the step size is reached";
    case {get_exitflag("MAXFUN_REACHED")}
        output.message = "The maximum number of function evaluations is reached";
    case {get_exitflag("FTARGET_REACHED")}
        output.message = "The target on the objective function is reached";
    case {get_exitflag("MAXIT_REACHED")}
        output.message = "The maximum number of iterations is reached";
    otherwise
        output.message = "Unknown exitflag";
end
