function [xval, fval, exitflag, output] = bds_powell(fun, x0, options)
%BDS Unconstrained nonlinear minimization (direct search with blocks).
%
%   XVAL = BDS(FUN, X0) starts at X0 and attempts to find
%   local minimizer X of the function FUN.  FUN is a function handle. FUN
%   accepts input X and returns a scalar function value F evaluated at X.
%   X0 should be a vector.
%
%   XVAL = BDS(FUN, X0, OPTIONS) minimizes with the
%   default optimization parameters replaced by values in the structure OPTIONS,
%   BLOCKWISE_DIRECT_SEARCH uses these options: num_blocks, MaxFunctionEvaluations, MaxFunctionEvaluations_factor,
%   expand, shrink, sufficient decrease factor, StepTolerance, ftarget, polling_inner,
%   blocks_strategy, with_cycling_memory, cycling, accept_simple_decrease.
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
%   num_blocks - number of blocks
%   MaxFunctionEvaluations - maximum of function evaluation
%   MaxFunctionEvaluations_factor - factor of maximum of function evaluation regarding to
%               dimensions.
%   expand - expanding factor of step size
%   shrink - shrinking factor of step size
%   sufficient_decrease_factor - factor of sufficient decrease condition
%   StepTolerance - StepTolerance of step size. If step size is below StepTolerance, then the
%        algorithm terminates.
%   ftarget - If function value is below ftarget, then the algorithm terminates.
%   polling_inner - polling strategy of indices in one block
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
    [fun] = verify_preconditions(fun, x0, options);
end

% The exit flag will be set at each possible exit of the algorithm.
% Therefore, if it is set to NaN on exit, it means that there is a bug.
exitflag = NaN;

% Transpose x0 if it is a row.
x0 = double(x0(:));

% Set the polling directions in D.
n = length(x0);
if strcmpi(options.Algorithm, "cbds") || strcmpi(options.Algorithm, "pbds")...
        || strcmpi(options.Algorithm, "ds") || strcmpi(options.Algorithm, "rbds")
    D = get_searching_set(n, options);
elseif strcmpi(options.Algorithm, "dspd")
    % Generate a vector which follows uniform distribution on the sphere of a unit ball.
    rv = NaN(n, 1);
    for i = 1:n
        rv(i) = randn(1);
    end
    Q = qr(rv);
    D = [Q, -Q];
end

% number of directions
num_directions = size(D, 2); 
% Set the default number of blocks.
if isfield(options, "num_blocks")
    num_blocks = options.num_blocks;
elseif strcmpi(options.Algorithm, "cbds") || strcmpi(options.Algorithm, "pbds")...
        || strcmpi(options.Algorithm, "rbds")
    % Default value is set as n, which is good for canonical with 2n directions. For
    % other situations, other value may be good.
    num_blocks = n;
elseif strcmpi(options.Algorithm, "dspd") || strcmpi(options.Algorithm, "ds")
    num_blocks = 1;
end

% If number of directions is less than number of blocks, then the number of
% blocks is defined as the number of directions.
num_blocks = min(num_directions, num_blocks);
% Default indices of blocks are 1:num_blocks.
block_indices = 1:num_blocks;

% Set MaxFunctionEvaluations to the maximum number of function evaluations. The default
% value is 1e5.
if isfield(options, "MaxFunctionEvaluations_factor") && isfield(options, "MaxFunctionEvaluations")
    MaxFunctionEvaluations = min(options.MaxFunctionEvaluations_factor*n, options.MaxFunctionEvaluations);
elseif isfield(options, "MaxFunctionEvaluations_factor")
    MaxFunctionEvaluations = options.MaxFunctionEvaluations_factor*n;
elseif isfield(options, "MaxFunctionEvaluations")
    MaxFunctionEvaluations = options.MaxFunctionEvaluations;
else
    MaxFunctionEvaluations = min(get_default_constant("MaxFunctionEvaluations"), get_default_constant("MaxFunctionEvaluations_factor")*n);
end

% Set the maximum of iterations. If complete polling is used, then the
% maximum number of iterations given below CANNOT be reached. If the
% opportunistic case is used, then the maximum number of iterations may
% be reached (although, we hope that it does not).
% ceil(10*MaxFunctionEvaluations/num_directions) may not be enough. Since there are some cases that
% maxit is exhausted and other terminations are not reached.
maxit = MaxFunctionEvaluations;

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

% Set default value of shuffling_period. Default value of shuffling_period
% should be set 1 since the algorithm visits all blocks for every iteration.
if strcmpi(options.Algorithm, "pbds") && isfield(options, "shuffling_period")
    shuffling_period = options.shuffling_period;
else
    shuffling_period = get_default_constant("shuffling_period");
end

% Set default value of replacement_delay. Default value of
% replacement_delay is set to be 0. Also, the value of replacement_delay
% should be less than or equal to num_blocks-1, otherwise, there will not exist
% such block going to be visited after num_blocks blocks have been visited.
if strcmpi(options.Algorithm, "rbds") && isfield(options, "replacement_delay")
    replacement_delay = min(options.replacement_delay, num_blocks-1);
else
    replacement_delay = min(get_default_constant("replacement_delay"), num_blocks-1);
end

% Set the default value for the boolean indicating whether the cycling
% strategy employed in the opportunistic case memorizes the history or not.
if isfield(options, "with_cycling_memory")
    with_cycling_memory = options.with_cycling_memory;
else
    with_cycling_memory = get_default_constant("with_cycling_memory");
end

% Set initial step size and alpha_hist to store the history of step size.
alpha_hist = NaN(num_blocks, maxit);
if isfield(options, "alpha_init")
    alpha_all = options.alpha_init*ones(num_blocks, 1);
else
    alpha_all = ones(num_blocks, 1);
end

% Divide the indices of the polling directions for each block.
searching_set_indices = divide_searching_set(num_directions, num_blocks);

% Initialize the computations.
% history of function values
fhist = NaN(1, MaxFunctionEvaluations);
% history of points having been visited
xhist = NaN(n, MaxFunctionEvaluations); 
% history of blocks having been visited
block_hist = NaN(1, MaxFunctionEvaluations);
xval = x0; % current iterate
fval = eval_fun(fun, xval);
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

if ~isfield(options, "powell_factor")
    powell_factor = get_default_constant("powell_factor");
else 
    powell_factor = options.powell_factor;
end

if ~isfield(options, "powell_factor_period")
    powell_factor_period = get_default_constant("powell_factor_period");
else 
    powell_factor_period = options.powell_factor_period;
end

% Start the actual computations.
% num_blocks blocks have been explored after the number of iteration goes from k to k+1.
for iter = 1 : maxit
    % record the value of alpha_all of the current iteration in alpha_hist.
    alpha_hist(:, iter) = alpha_all;
    
    % Why iter-1? Because the initial value of iter is 1 and when iter
    % increases by 1, the algorithm will visit num_blocks blocks when 
    % options.Algorithm = "pbds".
    if strcmpi(options.Algorithm, "pbds") && mod(iter - 1, shuffling_period) == 0
        % Make sure that `shuffling_period` is defined when `Algorithm` is "pbds".
        block_indices = randperm(num_blocks);
    end
    
    % Get the block that are going to be visited in this iteration.
    if strcmpi(options.Algorithm, "rbds")
        if replacement_delay == 0 || sum(~isnan(block_hist)) == 0
            block_indices = randi([1, num_blocks]);
        else
            % Recorder the number of blocks having been visited
            num_visited = sum(~isnan(block_hist));
            block_visited_slices_length = min(num_visited, replacement_delay);
            block_visited_slices = block_hist(num_visited-block_visited_slices_length+1:num_visited);
            % Set default value of block_indices
            block_initial_indices = 1:num_blocks;
            % Remove elements of block_indices appearing in block_visited_slice 
            block_real_indices = block_initial_indices(~ismember(block_initial_indices, block_visited_slices));
            % Produce a random index from block_real_indices.
            idx = randi(length(block_real_indices));
            block_indices = block_real_indices(idx);
        end
    end
    
    for i = 1:length(block_indices)
        % In case of permutation.
        i_real = block_indices(i);
        
        % Recorder the number of blocks having been visited
        num_visited = sum(~isnan(block_hist));
        % Update the block that going to be visited 
        block_hist(num_visited+1) = i_real;
        
        % get indices in the i-th block
        direction_indices = searching_set_indices{i_real}; 
        
        suboptions.MaxFunctionEvaluations = MaxFunctionEvaluations - nf;
        % Memory and cycling are needed since we permutate indices in inner_direct_search
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

    % Update alpha using Cunxin"s technique.
    if (max(alpha_all) <= powell_factor(1))
        powell_factor = shrink^powell_factor_period*powell_factor;
    else
        for i = 1:num_blocks
            if alpha_all(i)<=powell_factor(2)
                alpha_all(i) = powell_factor(2);
            end

        end
    end
    
end

% Set useful pieces on information about the solver"s history in output.
output.funcCount = nf;
output.fhist = fhist(1:nf);
output.xhist = xhist(:, 1:nf);
output.alpha_hist = alpha_hist(:, 1:min(iter, maxit));
% Recorder the number of blocks visited
num_blocks_visited = sum(~isnan(block_hist));
% Update the block that going to be visited
output.blocks_hist = block_hist(1:num_blocks_visited);


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

% Verify_postconditions: If debug_flag is true, then verify_postconditions is
% to verify output. If output_correctness is false, then assert will let code crash.
if debug_flag
    verify_postconditions(fun, xval, fval, exitflag, output);
end
