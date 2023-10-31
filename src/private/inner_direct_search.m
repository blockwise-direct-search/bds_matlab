function [xopt, fopt, exitflag, output] = inner_direct_search(fun, ...
    xopt, fopt, D, direction_indices, alpha, options)
%INNER_DIRECT_SEARCH performs a single iteration of classical direct search 
%   within a given block.
%
%   xopt = INNER_DIRECT_SEARCH(FUN, xopt, fopt, D, DIRECTION_INDICES, ALPHA, OPTIONS)
%   returns an xopt.
%
%   xopt = INNER_DIRECT_SEARCH(FUN, xopt, fopt, D, ...
%   DIRECTION_INDICES, ALPHA, OPTIONS) works with the structure OPTIONS, which includes
%   reduction_factor, ftarget, polling, with_cycling_memory, cycling.
%
%   [xopt, fopt] = INNER_DIRECT_SEARCH(...) returns the value of the objective function fopt 
%   at xopt.
%
%   [xopt, fopt, EXITFLAG] = INNER_DIRECT_SEARCH(...) returns an EXITFLAG that describes 
%   the exit condition.
%
%   [xopt, fopt, EXITFLAG, OUTPUT] = INNER_DIRECT_SEARCH(...) returns a structure OUTPUT including
%   funcCount, fhist, xhist, success, terminate, and direction_indices.
%
%   TERMINATE is initialized to be false. If terminate is updated to be true,
%   it means that either the number of function evaluations reaches maxfun, or ftarget is reached.
%
%   DIRECTION_INDICES is the indices of directions of this block in D.
%

% Set the value of sufficient decrease factor.
if ~isfield(options, "reduction_factor")
    reduction_factor = get_default_constant("reduction_factor");
else
    reduction_factor = options.reduction_factor;
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

% Explain why NaN is good. It is possible that this function returns
% with exitflag=NaN and this is NOT a bug. This is because other situations
% correspond to other normal values. Easy to see whether there is some bug 
% related to exitflag.
exitflag = NaN;

% Initialize some parameters before entering the loop.
n = length(xopt);
num_directions = length(direction_indices);
fhist = NaN(1, num_directions);
xhist = NaN(n, num_directions);
nf = 0; 
fbase = fopt;
xbase = xopt;
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
        xopt = xnew;
        fopt = fnew;
        terminate = true;
        information = "FTARGET_REACHED";
        exitflag = get_exitflag(information);
        break;
    end
    
    % Check whether the sufficient decrease condition is achieved.
    sufficient_decrease = (fnew + reduction_factor(3) * forcing_function(alpha)/2 < fbase);

    % Update the best point and the best function value.
    if fnew < fopt
        xopt = xnew;
        fopt = fnew;
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
output.direction_indices = direction_indices;
output.terminate = terminate;

end





