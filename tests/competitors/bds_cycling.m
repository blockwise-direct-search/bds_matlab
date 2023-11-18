function [xval, fval, exitflag, output] = bds_cycling(fun, x0, options)
%   bds_polling is designed under the framework of polling and cycling.
%
%   BLOCKWISE_DIRECT_SEARCH Unconstrained nonlinear minimization (direct search with blocks).
%
%   XVAL = BLOCKWISE_DIRECT_SEARCH(FUN, X0) starts at X0 and attempts to find
%   local minimizer X of the function FUN.  FUN is a function handle. FUN
%   accepts input X and returns a scalar function value F evaluated at X.
%   X0 should be a vector.
%
%   XVAL = BLOCKWISE_DIRECT_SEARCH(FUN, X0, OPTIONS) minimizes with the
%   default optimization parameters replaced by values in the structure OPTIONS,
%   BLOCKWISE_DIRECT_SEARCH uses these options: num_blocks, maxfun, maxfun_factor,
%   expand, shrink, sufficient decrease factor, StepTolerance, ftarget, polling_inner,
%   polling_outer, with_cycling_memory, cycling.
%
%   num_blocks - number of blocks
%   maxfun - maximum of function evaluation
%   maxfun_factor - factor of maximum of function evaluation regarding to
%               dimensions.
%   expand - expanding factor of step size
%   shrink - shrinking factor of step size
%   sufficient_decrease_factor - factor of sufficient decrease condition
%   StepTolerance - StepTolerance of step size. If step size is below StepTolerance, then the
%        algorithm terminates.
%   ftarget - If function value is below ftarget, then the algorithm terminates.
%   polling_inner - polling strategy of indices in one block
%   polling_outer - polling strategy of block_indices
%   with_cycling_memory - If with_cycling_memory is true, permutation will be executed on the array of last
%           iteration, otherwise, permutation will be executed on the initial array.
%   cycling - Possible values of cycling and the corresponding conditions
%   are listed below.
%
%   If index = 0, then there is no permutation.
%
%   0  No permutation.
%
%   1  The element of the index will be moved to the first element of array.
%
%   EXAMPLE
%   When array is 3 1 2 4 5, if index = 3, array will be 2 3 1 4 5 after 
%   cycling when with_cycling_memory is true; index will be 2, sort(index)
%   is 1 2 3 4 5 and array will be 2 1 3 4 5 after cycling when with_cycling_memory is false.
%
%   2  The element of the index and the following ones until end will be
%      moved ahead of array.
%
%   EXAMPLE
%   When array is 2 1 4 5 3, if index = 3, array will be 4 5 3 2 1 after cycling when 
%   with_cycling_memory is true; index will be 4, sort(index) is 1 2 3 4 5 and array will 
%   be 4 5 1 2 3 after cycling when with_cycling_memory is false.
%
%   3  The element of the following ones after index until end will be
%      moved ahead of array.
%
%   EXAMPLE
%   When array is 2 1 4 5 3 and index = 3, array will be 5 3 2 1 4 after cycling when 
%   with_cycling_memory is true; index will be 4, sort(index) is 1 2 3 4 5 and array 
%   will be 5 1 2 3 4 after cycling when with_cycling_memory is false.
%
%   4  The element of the following one after index will be moved ahead of array.
%
%   EXAMPLE
%   array is 4 1 2 3 5, if index = 3, array will be 3 4 1 2 5 after cycling when 
%   with_cycling_memory is true; index will be 2, sort(index) is 1 2 3 4 5 and array 
%   will be 3 1 2 4 5 after cycling when with_cycling_memory is false.
%
% [XVAL, FVAL] = BLOCKWISE_DIRECT_SEARCH(...) returns the value of the
% objective function, described in FUN, at XVAL.
%
% [XVAL, FVAL, EXITFLAG] = BLOCKWISE_DIRECT_SEARCH(...) returns an EXITFLAG
% that describes the exit condition. The information of EXITFLAG will be
% given in output.message.
%
% [XVAL, FVAL, EXITFLAG, OUTPUT] = BLOCKWISE_DIRECT_SEARCH(...) returns a
% structure OUTPUT with fields
%
%   fhist      History of function value
%   xhist      History of points that being calculated
%   funcCount  The number of function evaluations.
%
% the number of function evaluations in OUTPUT.funcCount,
% the history of function evaluation in OUTPUT.fhist, the history of points
% in OUTPUT.xhist and the l2-norm of gradient in OUTPUT.ghist (only for CUTEST
% PROBLEM).

% Set options to an empty structure if it is not supplied.
if nargin < 3
    options = struct();
end

% Precondition: If debug_flag is true, then pre-conditions is operated on
% input. If input_correctness is false, then assert may let the code crash.
debug_flag = is_debugging();
if debug_flag
    precondition_bds(fun, x0, options);
end

% The exit flag will be set at each possible exit of the algorithm.
% Therefore, if it is set to NaN on exit, it means that there is a bug.
exitflag = NaN;

% Transpose x0 if it is a row.
x0 = double(x0(:));

% Set the polling directions in D.
n = length(x0);
D = searching_set(n, options);
num_directions = size(D, 2); % number of directions

% Set the default number of blocks.
if isfield(options, "num_blocks")
    num_blocks = options.num_blocks;
else
    % TODO: this default value is good for canonical with 2n directions. For
    % other situations, other value may be good.
    num_blocks = n;
end

% If number of directions is less than number of blocks, then the number of
% blocks is defined as the number of directions.
num_blocks = min(num_directions, num_blocks);

% Set maxfun to the maximum number of function evaluations. The default
% value is 1e4.
if isfield(options, "maxfun_factor")
    maxfun = options.maxfun_factor*n;
    if isfield(options, "maxfun")
        maxfun = min(options.maxfun, maxfun);
    end
elseif isfield(options, "maxfun")
    maxfun = options.maxfun;
else
    maxfun = get_default_constant("maxfun");
end

% Set the maximum of iterations. If complete polling is used, then the
% maximum number of iterations given below CANNOT be reached. If the
% opportunistic case is used, then the maximum number of iterations may
% be reached (although, we hope that it does not).
% ceil(10*maxfun/num_directions) may not be enough. Since there are some cases that
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
    % TODO: Check whether a "smarter" value is not possible, such as
    % "10 * eps * n" for example.
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

% Set the default outer polling strategy. This is the polling strategy
% employed to cycle the blocks.
if ~isfield(options, "polling_outer")
    options.polling_outer = get_default_constant("polling_outer");
end

% Set the default cycling strategy. The variable CANNOT be named "cycling,"
% because there is a function named "cycling."
if isfield(options, "cycling_outer")
    cycling_outer = options.cycling_outer;
else
    cycling_outer = get_default_constant("cycling_outer");
end

if isfield(options, "cycling_inner")
    cycling_inner = options.cycling_inner;
else
    cycling_inner = get_default_constant("cycling_inner");
end


% Set the default value for the boolean indicating whether the cycling
% strategy employed in the opportunistic case memorizes the history or not.
if isfield(options, "with_cycling_memory")
    with_cycling_memory = options.with_cycling_memory;
else
    with_cycling_memory = get_default_constant("with_cycling_memory");
end

% Set the initial block indices and the corresponding initial step sizes.
block_indices = 1:num_blocks;

if isfield(options, "alpha_init")
    alpha_all = options.alpha_init*ones(num_blocks, 1);
else
    alpha_all = ones(num_blocks, 1);
end

% Divide the indices of the polling directions for each block.
% TODO: Tell Zaikun that Tom disagrees with this name.
searching_set_indices = divide_searching_set(num_directions, num_blocks);

% Initialize the computations.
fhist = NaN(1, maxfun); % history of function values
xhist = NaN(n, maxfun); % history of iterates
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

% Initialize success_block_index. It is the index of the block where a
% success occurred and is hence a value between 1 and num_blocks. When no
% sufficient decrease is observed for all directions in all blocks, we set
% success_block_index to -1, which will be more easier to port python and C in
% the future.
% success_block_index = -1;

% Start the actual computations.
for k = 1 : maxit

    % Let xbase be the point from which the polling directions are
    % employed. In one iteration, all the block use the same base point.
    % The corresponding value of the objective function is stored in fbase.
    xbase = xval(:);
    fbase = fval;

    % We need to ensure that success_index = -1 before entering the next
    % iteration, otherwise success_index will be inherited from last iteration, which is wrong.
    success_block_index = -1;

    terminate = false;
    % Do the computations among the blocks.
    for i = 1 : length(block_indices)
        i_real = block_indices(i);
        direction_indices = searching_set_indices{i_real}; % acquire indices in block i_real

        suboptions.maxfun = maxfun - nf;
        % Memory and cycling are needed since we permutate indices in inner_direct_search
        suboptions.cycling = cycling_inner;
        suboptions.with_cycling_memory = with_cycling_memory;
        suboptions.sufficient_decrease_factor = sufficient_decrease_factor;
        suboptions.ftarget = ftarget;
        suboptions.polling_inner = options.polling_inner;
        suboptions.accept_simple_decrease = accept_simple_decrease;

        [xval, fval, sub_exitflag, suboutput] = inner_direct_search(fun, xval,...
            fval, xbase, fbase, D(:, direction_indices), direction_indices,...
            alpha_all(i_real), suboptions);

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

        % Update the step sizes.
        if success
            alpha_all(i_real) = expand * alpha_all(i_real);
        else
            alpha_all(i_real) = shrink * alpha_all(i_real);
        end

        % Terminate the computations if the largest step size is below a
        % given StepTolerance.
        % TODO: Is it normal to check whether "SMALL_ALPHA" is reached
        % directly after updating the step sizes, or should be do one more
        % iteration with the last value of the step sizes?
        if max(alpha_all) < alpha_tol
            terminate = true;
            exitflag = get_exitflag("SMALL_ALPHA");
            break
        end

        % Store the index of the block in which a direction provided a
        % sufficient decrease (if any). This direction will be used to
        % cycle the block indices.
        if success && ~strcmpi(options.polling_outer, "complete")
            success_block_index = i;
            break;
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
    if k == maxit
        exitflag = get_exitflag("MAXIT_REACHED");
    end

    % Cycle the block indices in the opportunistic case, following the
    % strategy given in options.polling_outer, cycling, and with_cycling_memory.
    if ~strcmpi(options.polling_outer, "complete")
        block_indices = cycling(block_indices, success_block_index, cycling_outer, with_cycling_memory);
    end

end

% Set useful pieces on information about the solver"s history in output.
output.funcCount = nf;
output.fhist = fhist(1:nf);
output.xhist = xhist(:, 1:nf);
%output.alpha_hist = alpha_hist;

% Postcondition: If debug_flag is true, then post-conditions is operated on
% output. If output_correctness is false, then assert will let code crash.
if debug_flag
    postcondition_bds(fun, xval, fval, exitflag, output);
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
