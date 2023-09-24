function [xval, fval, exitflag, output] = dspd(fun, x0, options)
%DSPD (direct search probabilistic descent) solves unconstrained optimization
%   problems without using derivatives. 
%
%   TODO: here(MATLAB version for dspd)
%
%   XVAL = DSPD(FUN, X0) returns an approximate minimizer XVAL of the function handle FUN, starting the
%   calculations at X0. FUN must accept input X and returns a scalar, which is the function value
%   evaluated at X. X0 should be a vector.
%
%   XVAL = DSPD(FUN, X0, OPTIONS) performs the computations with the options in OPTIONS. It should be a
%   structure, with the following fields:
%   
%   maxfun                      Maximum of function evaluations.
%   maxfun_factor                  Factor to define maximum number of function evaluations as a multiplier
%                               of the dimension of the problem.    
%   expand                      Expanding factor of step size.
%   shrink                      Shrinking factor of step size.
%   sufficient_decrease_factor  Factor of sufficient decrease condition.
%   StepTolerance               The tolerance for testing whether the step size is small enough.
%   ftarget                     Target of function value. If function value is below ftarget, 
%                               then the algorithm terminates.
%   polling_inner               Polling strategy of each block.
%   with_cycling_memory         In the opportunistic case (polling_inner == "opportunistic"), 
%                               with_memory decides whether the cycling strategy memorizes 
%                               the history or not.
%   cycling_inner               Cycling strategy employed in the opportunistic case.
%   accept_simple_decrease      Whether the algorithm accepts simple decrease or not. 
%   seed                        Only used by randomized strategy for reproducibility.
%
%   [XVAL, FVAL] = DSPD(...) also returns the value of the objective function FUN at the 
%   solution XVAL.
%
%   [XVAL, FVAL, EXITFLAG] = DSPD(...) returns an EXITFLAG that describes the exit 
%   condition. The possible values of EXITFLAG are 0, 1, 2, 3.
%
%   0    The StepTolerance of the step size is reached.
%   1    The maximum number of function evaluations is reached.
%   2    The target of the objective function is reached.
%   3    The maximum number of iterations is reached.
%
%   [XVAL, FVAL, EXITFLAG, OUTPUT] = DSPD(...) returns a
%   structure OUTPUT with the following fields: 
%
%   fhist        History of function values.
%   xhist        History of points visited.
%   alpha_hist   History of step size.
%   funcCount    The number of function evaluations.
%   message      The information of EXITFLAG.
%
%   Copyright 2023, Haitian Li and Zaikun Zhang.
%   All rights reserved.
%

% Get the information of current MATLAB.
matlabVersion = ver('MATLAB');

% Extract the main version and the minor version of MATLAB.
versionParts = strsplit(matlabVersion.Version, '.');

% Change the version of matlab into floating number.
majorVersion = str2double(versionParts{1});
minorVersion = str2double(versionParts{2});

% Check whether the version of MATLAB is released lower than 2017a.
if majorVersion < 9 || (majorVersion == 9 && minorVersion < 1)
    error("The version of MATLAB is low, please update to 2017a or later versions")
end

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

% Set the value of expand factor.
if isfield(options, "expand")
    expand = options.expand;
else
    expand = get_default_constant("expand");
end

% Under dspd, expand factor could not be 1.
if expand == 1
    error("If expanding factor is set to be 1, convergence may not be ensured")
end

% Set the value of shrink factor.
if isfield(options, "shrink")
    shrink = options.shrink;
else
    shrink = get_default_constant("shrink");
end

% Get the number of directions.
if isfield(options, "num_random_vectors")
    m = max(options.num_random_vectors, ceil(log2(1-log(shrink))/log(expand)));
else
    m = max(get_default_constant("num_random_vectors"), ceil(log2(1-log(shrink))/log(expand)));
end

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

% Set the value of sufficient decrease factor.
if isfield(options, "sufficient_decrease_factor")
    sufficient_decrease_factor = options.sufficient_decrease_factor;
else
    sufficient_decrease_factor = get_default_constant("sufficient_decrease_factor");
end

% Set the boolean value of accept_simple_decrease, which is for updating xval and fval, but not 
% for stepsize. 
if isfield(options, "accept_simple_decrease")
    accept_simple_decrease = options.accept_simple_decrease;
else
    accept_simple_decrease = get_default_constant("accept_simple_decrease");
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
% TODO: write it more detailedly.
if isfield(options, "with_cycling_memory")
    with_cycling_memory = options.with_cycling_memory;
else
    with_cycling_memory = get_default_constant("with_cycling_memory");
end

% Initialize the step sizes and alpha_hist, which is the history of step sizes.
alpha_hist = NaN(1, maxit);
if isfield(options, "alpha_init")
    alpha = options.alpha_init;
else
    alpha = 1;
end

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
    % Record the value of alpha_all of the current iteration in alpha_hist.
    alpha_hist(iter) = alpha;

    % Generate the searching set whose directions are uniformly distributed on the unit sphere
    % for each iteration when options.Algorithm is "dspd".
    if m == 1
        rv = randn(n, 1);
        % Normalize rv.
        rv = rv ./ norm(rv);
        D = [rv, -rv];
    else
        D = randn(n, m);
        % Normalize D. vecnorm is introduced in MATLAB 2017a for the first time.
        % Here we are dividing a matrix by a row vector using implicit expansion.
        D = D ./ vecnorm(D);
    end

    % Get indices of directions in the i-th block.
    direction_indices = 1:size(D, 2);
   
    suboptions.maxfun = maxfun - nf;
    suboptions.cycling = cycling_inner;
    suboptions.with_cycling_memory = with_cycling_memory;
    suboptions.sufficient_decrease_factor = sufficient_decrease_factor;
    suboptions.ftarget = ftarget;
    suboptions.polling_inner = options.polling_inner;
    suboptions.accept_simple_decrease = accept_simple_decrease;

    [xval, fval, sub_exitflag, suboutput] = inner_direct_search(fun, xval,...
        fval, D, direction_indices, alpha, suboptions);

    % Update the history of step size.
    alpha_hist(iter) = alpha;

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

    % Check whether a sufficient decrease has been achieved in inner_direct_search.
    success = suboutput.success;

    % Update the step sizes and store the history of step sizes.
    if success
        alpha = expand * alpha;
    else
        alpha = shrink * alpha;
    end

    % Terminate the computations if the largest component of step size is below a
    % given StepTolerance.
    if alpha < alpha_tol
        exitflag = get_exitflag("SMALL_ALPHA");
        break
    end
end

% Check whether MAXIT is reached.
if iter == maxit
    exitflag = get_exitflag("MAXIT_REACHED");
end

% Truncate HISTORY into an nf length vector.
output.funcCount = nf;
output.fhist = fhist(1:nf);
output.xhist = xhist(:, 1:nf);
output.alpha_hist = alpha_hist(1:min(iter, maxit));
 
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

% verify_postconditions is to detect whether output is in right form when debug_flag is true.
if debug_flag
    verify_postconditions(fun, xval, fval, exitflag, output);
end

end
