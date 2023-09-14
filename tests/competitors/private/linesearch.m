function [xval, fval, exitflag, output] = linesearch(fun, ...
    xval, fval, D, direction_indices, alpha, options)
%LINESEARCH peforms a single iteration of classical direct search 
%   within a given block.
%
%   XVAL = LINESEARCH(FUN, XVAL, FVAL, D, DIRECTION_INDICES, ALPHA, OPTIONS)
%   returns a XVAL.
%
%   XVAL = LINESEARCH(FUN, XVAL, FVAL, D, ...
%   DIRECTION_INDICES, ALPHA, OPTIONS) works with the structure OPTIONS, which includes
%   sufficient_decrease_factor, ftarget.
%
%   [XVAL, FVAL] = LINESEARCH(...) returns the value of the objective function FVAL 
%   at XVAL.
%
%   [XVAL, FVAL, EXITFLAG] = LINESEARCH(...) returns an EXITFLAG that describes 
%   the exit condition.
%
%   [XVAL, FVAL, EXITFLAG, OUTPUT] = LINESEARCH(...) returns a structure OUTPUT including
%   funcCount, fhist, xhist, success, terminate.
%
%   SUCCESS is initialized to be false. If success is updated to be true,
%   it means that there at least exists some direction satisfying sufficient decrease.
%
%   TERMINATE is initialized to be false. If terminate is updated to be true,
%   it means that either the number of function evaluations reaches maxfun, or ftarget is reached.
%
%   DIRECTION_INDICES is indices of directions of this block in D.
%

% Set the value of sufficient decrease factor.
if ~isfield(options, "sufficient_decrease_factor")
    sufficient_decrease_factor = get_default_constant("sufficient_decrease_factor");
else
    sufficient_decrease_factor = options.sufficient_decrease_factor;
end

% Set the value of expanding factor.
expand = options.expand;

% Set the value of cycling_strategy, which represents the cycling strategy inside each block.
cycling_strategy = 1;

% Set the boolean value of WITH_CYCLING_MEMORY. 
if isfield(options, "with_cycling_memory")
    with_cycling_memory = options.with_cycling_memory;
else
    with_cycling_memory = get_default_constant("with_cycling_memory");
end

% Set ftarget of objective function.
if ~isfield(options, "ftarget")
    ftarget = get_default_constant("ftarget");    
else
    ftarget = options.ftarget;
end

% Explain why NaN is good. It is possible that this function returns
% with exitflag=NaN and this is NOT a bug. This is because other situatons
% are corresponding to other normal values. Easy to see whether there is
% some bug related to exitflag.
exitflag = NaN;

% Initialize some parameters before entering the loop.
n = length(xval);
num_directions = length(direction_indices);
fhist = NaN(1, options.maxfun);
xhist = NaN(n, options.maxfun);
success = false;
nf = 0; 
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
    sufficient_decrease = (fnew + sufficient_decrease_factor * alpha^2/2 < fbase);
    
    % if sufficient decrease
    if sufficient_decrease && fnew < fval
        fval = fnew;
        xval = xnew;
    end

    success = sufficient_decrease;

    while sufficient_decrease
        % error("This part of code is not used.")
        alpha = alpha*expand;
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

        % Stop the loop if no more function evaluations can be performed. 
        % Note that this should be checked before evaluating the objective function.
        if nf >= options.maxfun
            terminate = true;
            exitflag = get_exitflag("MAXFUN_REACHED");
            break;
        end
        
        if strcmpi(options.linesearch_type, "standard")
            sufficient_decrease = (fnew + sufficient_decrease_factor * alpha^2 < fbase);
        elseif strcmpi(options.linesearch_type, "new")
            sufficient_decrease = (fnew + sufficient_decrease_factor * ((expand-1)...
                *alpha)^2 < fhist(nf-1));    
        end

        

        if sufficient_decrease
            fval = fnew;
            xval = xnew;
        end
    end
     
    if success
        direction_indices = cycling(direction_indices, j, cycling_strategy, with_cycling_memory);
        break;
    end
end

% Truncate FHIST and XHIST into an nf length vector.
output.fhist = fhist(1:nf);
output.xhist = xhist(:, 1:nf);
output.nf = nf;
output.success = success;
output.direction_indices = direction_indices;
output.terminate = terminate;
output.stepsize = alpha;

end




