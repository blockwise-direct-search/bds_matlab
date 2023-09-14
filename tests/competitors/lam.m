function [xval, fval, exitflag, output] = lam(fun, x0, options)
%LAM (Linesearch Algorithm Model) solves unconstrained optimization problems without using derivatives. 
%
%
%   XVAL = LAM(FUN, X0) returns an approximate minimizer XVAL of the function handle FUN, starting the
%   calculations at X0. FUN must accept input X and returns a scalar, which is the function value
%   evaluated at X. X0 should be a vector.
%
%   XVAL = LAM(FUN, X0, OPTIONS) performs the computations with the options in OPTIONS. It should be a
%   structure, with the following fields:
%   
%   nb                          Number of blocks.
%   maxfun                      Maximum of function evaluations.
%   maxfun_dim                  Factor to define maximum number of function evaluations as a multiplier
%                               of the dimension of the problem.    
%   expand                      Expanding factor of step size (not the same with direct search).
%   shrink                      Shrinking factor of step size.
%   sufficient_decrease_factor  Factor of sufficient decrease condition.
%   StepTolerance               The tolerance for testing whether the step size is small enough.
%   ftarget                     Target of function value. If function value is below ftarget, 
%                               then the algorithm terminates.
%   stepsize_factor             constant for adjusting stepsize before computing.
%
%   [XVAL, FVAL] = LAM(...) also returns the value of the objective function FUN at the 
%   solution XVAL.
%
%   [XVAL, FVAL, EXITFLAG] = LAM(...) returns an EXITFLAG that describes the exit 
%   condition. The possible values of EXITFLAG are 0, 1, 2, 3.
%
%   0    The StepTolerance of the step size is reached.
%   1    The maximum number of function evaluations is reached.
%   2    The target of the objective function is reached.
%   3    The maximum number of iterations is reached.
%
%   [XVAL, FVAL, EXITFLAG, OUTPUT] = LAM(...) returns a
%   structure OUTPUT with the following fields: 
%
%   fhist        History of function values.
%   xhist        History of points visited.
%   funcCount    The number of function evaluations.
%   message      The information of EXITFLAG.
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

% Get the dimension of the problem.
n = length(x0);

% Get the searching set, the number of directions and blocks respectively.
%options.direction = "canonical";
D = get_searching_set(n, options);
m = size(D, 2);
nb = n;

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

% Set the value of StepTolerance. The algorithm will terminate if the stepsize is less than 
% the StepTolerance.
if isfield(options, "StepTolerance")
    alpha_tol = options.StepTolerance;
else
    alpha_tol = get_default_constant("StepTolerance");
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

% Set the boolean value of WITH_CYCLING_MEMORY. 
if isfield(options, "with_cycling_memory")
    with_cycling_memory = options.with_cycling_memory;
else
    with_cycling_memory = get_default_constant("with_cycling_memory");
end

% Set the value of stepsize_factor.
% if isfield(options, "stepsize_factor")
%     stepsize_factor = options.stepsize_factor;
% else
%     stepsize_factor = 0;
% end

% Set the target of the objective function.
if isfield(options, "ftarget")
    ftarget = options.ftarget;
else
    ftarget = get_default_constant("ftarget");
end

% Initialize the step sizes.
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

% Start the actual computations.
for iter = 1:maxit

    %alpha_max = max(alpha_all); 

    for i = 1:length(block_indices)
        % If block_indices is 1 3 2, then block_indices(2) = 3, which is the real block that we are
        % going to visit.
        i_real = block_indices(i);
        
        alpha = alpha_all(i_real);
        %alpha_bar = max(alpha_all(i_real), stepsize_factor*alpha_max);

        % Get indices of directions in the i-th block.
        direction_indices = searching_set_indices{i_real}; 
        
        suboptions.maxfun = maxfun - nf;
        suboptions.sufficient_decrease_factor = sufficient_decrease_factor;
        suboptions.with_cycling_memory = with_cycling_memory;
        %suboptions.expand = expand;
        suboptions.ftarget = ftarget;
        
        [xval, fval, sub_exitflag, suboutput] = linesearch(fun, xval,...
            fval, D(:, direction_indices), direction_indices,...
            alpha, suboptions);
        
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
        
        % Update the step sizes.
        if success
            % alpha_all(i_real) = suboutput.stepsize;
            % alpha_all(i_real) = expand * suboutput.stepsize;
            % alpha_all(i_real) = expand * alpha_bar;
            alpha_all(i_real) = expand * alpha;
        else
            % alpha_all(i_real) = shrink * alpha_bar;
            alpha_all(i_real) = shrink * alpha;
            %alpha_all(i_real) = shrink * alpha_all(i_real);
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

% Truncate HISTORY into an nf length vector.
output.funcCount = nf;
output.fhist = fhist(1:nf);
output.xhist = xhist(:, 1:nf);

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

