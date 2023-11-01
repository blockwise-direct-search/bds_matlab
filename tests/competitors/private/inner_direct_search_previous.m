function [xval, fval, exitflag, output] = inner_direct_search_previous(fun, ...
    xval, fval, D, direction_indices, alpha, options)
%INNER_DIRECT_SEARCH performs a single iteration of classical direct search
%   within a given block.
%
%   XVAL = INNER_DIRECT_SEARCH(FUN, XVAL, FVAL, D, DIRECTION_INDICES, ALPHA, OPTIONS)
%   returns an XVAL.
%
%   XVAL = INNER_DIRECT_SEARCH(FUN, XVAL, FVAL, D, ...
%   DIRECTION_INDICES, ALPHA, OPTIONS) works with the structure OPTIONS, which includes
%   sufficient_decrease_factor, ftarget, polling, with_cycling_memory, cycling.
%
%   [XVAL, FVAL] = INNER_DIRECT_SEARCH(...) returns the value of the objective function FVAL
%   at XVAL.
%
%   [XVAL, FVAL, EXITFLAG] = INNER_DIRECT_SEARCH(...) returns an EXITFLAG that describes
%   the exit condition.
%
%   [XVAL, FVAL, EXITFLAG, OUTPUT] = INNER_DIRECT_SEARCH(...) returns a structure OUTPUT including
%   funcCount, fhist, xhist, success, terminate, and direction_indices.
%
%   SUCCESS is initialized to be false. If success is updated to be true,
%   it means that there at least exists some direction satisfying sufficient decrease.
%
%   REDUCTION is initialize to be false. If reduction is updated to be true,
%   it means that there at least exists some direction achieving
%   reduction_ratio. In this case, we will not shrink the step size.
%
%   TERMINATE is initialized to be false. If terminate is updated to be true,
%   it means that either the number of function evaluations reaches maxfun, or ftarget is reached.
%
%   DIRECTION_INDICES is the indices of directions of this block in D.
%

% Set the value of sufficient decrease factor.
if ~isfield(options, "sufficient_decrease_factor")
    sufficient_decrease_factor = get_default_constant("sufficient_decrease_factor");
else
    sufficient_decrease_factor = options.sufficient_decrease_factor;
end

% Set target of objective function.
if ~isfield(options, "ftarget")
    ftarget = get_default_constant("ftarget");
else
    ftarget = options.ftarget;
end

% Set the value of polling_inner. This is the polling strategy employed within one block.
if isfield(options, "polling_inner")
    polling_inner = options.polling_inner;
else
    polling_inner = get_default_constant("polling_inner");
end

% Set the value of cycling_inner, which represents the cycling strategy inside each block.
if isfield(options, "cycling_inner")
    cycling_strategy = options.cycling_inner;
else
    cycling_strategy = get_default_constant("cycling_inner");
end

% Set the boolean value of WITH_CYCLING_MEMORY.
if isfield(options, "with_cycling_memory")
    with_cycling_memory = options.with_cycling_memory;
else
    with_cycling_memory = get_default_constant("with_cycling_memory");
end

% Set the forcing function, which is a function of alpha.
if isfield(options, "forcing_function")
    forcing_function = options.forcing_function;
else
    forcing_function = get_default_constant("forcing_function");
end

% Explain why NaN is good. It is possible that this function returns
% with exitflag=NaN and this is NOT a bug. This is because other situations
% correspond to other normal values. Easy to see whether there is some bug
% related to exitflag.
exitflag = NaN;

% Initialize some parameters before entering the loop.
n = length(xval);
num_directions = length(direction_indices);
fhist = NaN(1, num_directions);
xhist = NaN(n, num_directions);
nf = 0;
success = false;
reduction = false;
fbase = fval;
xbase = xval;
terminate = false;

for j = 1 : num_directions
    % Stop the loop if no more function evaluations can be performed.
    % Note that this should be checked before evaluating the objective function.
    if nf >= options.maxfun
        terminate = true;
        exitflag = get_exitflag("MAXFUN_REACHED");
        break;
    end

    % Evaluate the objective function for the current polling direction.
    xnew = xbase+alpha*D(:, j);
    fnew = eval_fun(fun, xnew);
    nf = nf+1;
    fhist(nf) = fnew;
    xhist(:, nf) = xnew;

    % Stop the computations once the target value of the objective function
    % is achieved.
    if fnew <= ftarget
        xval = xnew;
        fval = fnew;
        terminate = true;
        information = "FTARGET_REACHED";
        exitflag = get_exitflag(information);
        break;
    end

    % Check whether the sufficient decrease condition is achieved.
    if sufficient_decrease_factor <= 0
        sufficient_decrease = (fnew < fbase);
    else
        sufficient_decrease = (fnew + sufficient_decrease_factor(3) * forcing_function(alpha) < fbase);
        if ~sufficient_decrease
            % Check whether the smaller reduction_ratio is achieved. If
            % reduction is set to be false always, then the framework
            % is the same as the framework with only one
            % sufficient_decrease_factor. An effective way is to set
            % sufficient_decrease_factor(1) large enough to make reduction
            % always be false. Besides that, sufficient_decrease_factor(3) should be greater than
            % or equal to sufficient_decrease_factor(1) to make sure that the sufficient decrease
            % condition above is stronger than the sufficient decrease condition below.
            % reduction = false.
            reduction = (reduction || fnew + sufficient_decrease_factor(2) * forcing_function(alpha) < fbase);
        end
    end

    % Success is initialized to be false. Once there exists some direction satisfying sufficient
    % decrease, success will always be true inside this function.
    success = (success || sufficient_decrease);

    % We will update the value of xval and fval if the sufficient decrease condition below is achieved
    % and fnew < fval. This is because we want to make sure that the value of xval and fval is always
    % updated to be the best value. When we set sufficient_decrease_factor(1) to be 0, this is the same
    % as accept simple decrease. Also, sufficient_decrease_factor(3) should be greater than sufficient_decrease_factor(1)
    % to make sure that the sufficient decrease condition below is weaker than the sufficient decrease condition above.
    if (fnew + sufficient_decrease_factor(1) * alpha^3/2 < fbase) && fnew < fval
        xval = xnew;
        fval = fnew;
    end

    % In the opportunistic case, if the current iteration achieves sufficient decrease,
    % stop the computations after cycling the indices of the polling directions. The reason
    % why we cycle indices here is because inner_direct_search is called in a loop in outer_direct_search.
    if sufficient_decrease && ~strcmpi(polling_inner, "complete")
        direction_indices = cycling(direction_indices, j, cycling_strategy, with_cycling_memory);
        break;
    end
end

% Truncate FHIST and XHIST into a vector of length nf.
output.fhist = fhist(1:nf);
output.xhist = xhist(:, 1:nf);
output.nf = nf;
output.success = success;
output.reduction = reduction;
output.direction_indices = direction_indices;
output.terminate = terminate;

end



