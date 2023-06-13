function [xval, fval, exitflag, output] = direct_search(fun, x0, options)

% direct_search Unconstrained nonlinear minimization (direct search without blocks).

% XVAL = DIRECT_SEARCH(FUN, X0) starts at X0 and attempts to find local
% minimizer X of the function FUN.  FUN is a function handle. FUN
% accepts input X and returns a scalar function value F evaluated at X.
% X0 should be a vector.
%
% XVAL = DIRECT_SEARCH(FUN, X0, OPTIONS) minimizes with the default
% optimization parameters replaced by values in the structure OPTIONS.
% DIRECT_SEARCH uses these options: nb, maxfun, maxfun_dim,
% expand, shrink, sufficient decrease factor, StepTolerance, ftarget, polling,
% with_memory, cycling.
%
% [XVAL, FVAL] = DIRECT_SEARCH(...) returns the value of the objective
% function, described in FUN, at XVAL.
%
% [XVAL, FVAL, EXITFLAG] = DIRECT_SEARCH(...) returns an EXITFLAG
% that describes the exit condition.
%
% [XVAL, FVAL, EXITFLAG, OUTPUT] = DIRECT_SEARCH(...) returns a OUTPUT
% with the number of function evaluations in OUTPUT.funcCount, history of
% function evaluation in OUTPUT.fhist, the history of points in
% OUTPUT.xhist and the l2-norm of gradient in OUTPUT.ghist (only for CUTEST
% PROBLEM).

% Set options to an empty structure if it is not supplied.
if nargin < 3
    options = struct();
end

% The exit flag will be set at each possible exit of the algorithm.
% Therefore, if it is set to NaN on exit, it means that there is a bug.
exitflag = NaN;

% The initial point should be a column vector (or even better, transpose
% it if it is not the case).
num_column = size(x0,2);
if num_column ~= 1
    x0 = x0';
end

% Set the polling directions in D.
n = length(x0);
D = searching_set(n, options);
m = size(D, 2); % number of directions

% Set maxfun to the maximum number of function evaluations. The default
% value is 1e4.
if isfield(options, "maxfun_dim")
    maxfun = options.maxfun_dim*n;
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
maxit = ceil(10*maxfun/m);

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

% Set the default polling strategy. This is the polling strategy
% employed within one block.
if ~isfield(options, "polling")
    options.polling = options.polling_inner;
% else
%     options.polling = get_default_constant("polling");
end

% Set the default cycling strategy. The variable CANNOT be named "cycling,"
% because there is a function named "cycling."
if isfield(options, "cycling")
    cycling_strategy = options.cycling;
else
    cycling_strategy = get_default_constant("cycling");
end

% Set the default value for the boolean indicating whether the cycling
% strategy employed in the opportunistic case memorizes the history or not.
if isfield(options, "with_memory")
    with_memory = options.with_memory;
else
    with_memory = get_default_constant("with_memory");
end

% Set the initial indices and initial step sizes.
indices = 1:m;

% Set initial step size
if isfield(options, "alpha0")
    alpha = options.alpha0;
else
    alpha = 1;
end

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
if  fval < ftarget
    terminate = true;
    information = "FTARGET_REACHED";
    exitflag = get_exitflag(information);
else
    terminate = false;
end

% Initialize success. It is the index of the indices where a
% success occurred and is hence a value between 1 and m. When no
% sufficient decrease is observed for all directions, we set
% success_index to -1, which will be more easier to port into python or C.
success_index = -1;
for k = 1 : maxit

    % If fval is small than ftarget by using only one function evaluation, then
    % there is no need to compute any more. Except above situation, the
    % code paragraph below will never execute break.
    if terminate
        break;
    end

    % Initialize success and only update when a direction provided a sufficient
    % decrease with opportunistic polling.
    success = false;

    % Let xbase be the point from which the polling directions are
    % employed. In one iteration, all the block use the same base point.
    % The corresponding value of the objective function is stored in fbase.
    xbase = xval;
    fbase = fval;

    % Cycle the indices in the opportunistic case, following the
    % strategy given in options.polling, cycling, and with_memory.
    if ~strcmpi(options.polling, "complete")
        indices = cycling(indices, success_index, cycling_strategy, with_memory);
    end

    % We need to ensure that success_index = 0 before entering the next
    % iteration, otherwise success_index will be inherited from last iteration, which is wrong.
    success_index = -1;
    for i = 1 : m
        i_real = indices(i); % acquire indices
        xnew = xbase+alpha*D(:,i_real);
        fnew = fun(xnew);
        nf = nf + 1;
        fhist(nf) = fnew;
        xhist(:,nf) = xnew;

        % FIXME: Later on, once every work, the following "if" can be
        % placed after the current "if", so that every simple decrease
        % is accepted (without requiring sufficient decrease).
        % If we want to accept simple decrease, just use the condition of fnew < fval to test.
        % TODO:
        % 1. Comment why the following line is wrong.
        % 2. What if we update fnew and xnew whenever there is a smple decrease?
        %success = (fnew <= fbase - sufficient_decrease_factor * alpha^2 / 2);
        if fnew <= fbase - sufficient_decrease_factor * alpha^2 / 2
            success = true;
            if fnew < fval
                xval = xnew;
                fval = fnew;
            end
        end

        % If the maximum number of function evaluations is reached, then
        % the current loop (among the searching directions) is broken.
        % Then, alpha will be updated anyway (for simplicity of the code),
        % and the outer loop will be broken.
        if nf >= maxfun
            information = "MAXFUN_REACHED";
            terminate = true;
            exitflag = get_exitflag(information);
            break;
        end

        if fval <= ftarget
            information = "FTARGET_REACHED";
            terminate = true;
            exitflag = get_exitflag(information);
            break;
        end

        % If this section is reached, then both the maximum number of
        % function evaluations and the target value on the objective
        % function have not been reached. We then entertain a normal
        % iteration.
        if success && ~strcmpi(options.polling, "complete")
            success_index = i;
            break;
        end

    end

    % Update the step size.
    if success
        alpha = expand * alpha;
    else
        alpha = shrink * alpha;
    end

    % Terminate the computations if the largest step size is below a
    % given StepTolerance.
    if alpha < alpha_tol
        information = "SMALL_ALPHA";
        terminate = true;
        exitflag = get_exitflag(information);
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

end

% Define the output
output.funcCount = nf;
output.fhist = fhist(1:nf);
output.xhist = xhist(:, 1:nf);
switch exitflag
    case {get_exitflag("SMALL_ALPHA")}
        output.message = "The StepTolerance on the step size is reached";
    case {get_exitflag("MAXFUN_REACHED")}
        output.message = "Maximum number of function evaluation performed";
    case {get_exitflag("FTARGET_REACHED")}
        output.message = "Target on the objective function is reached";
    case {get_exitflag("MAXIT_REACHED")}
        output.message = "Maximum number of iterations performed";
    otherwise
        output.message = "Unknown exitflag";
end
