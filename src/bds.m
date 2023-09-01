function [xval, fval, exitflag, output] = bds(fun, x0, options)
%BDS (blockwise direct search) solves unconstrained optimization problems without using derivatives. 
%
%   XVAL = BDS(FUN, X0) starts at X0 and attempts to find
%   local minimizer X of the function FUN.  FUN is a function handle. FUN
%   accepts input X and returns a scalar function value F evaluated at X.
%   X0 should be a vector.
%
%   XVAL = BDS(FUN, X0, OPTIONS) minimizes with the
%   default optimization parameters replaced by values in the structure OPTIONS.
%   OPTIONS includes nb, maxfun, maxfun_dim, expand, shrink, sufficient decrease factor, 
%   StepTolerance, ftarget, polling_inner, with_memory, cycling, 
%   accept_simple_decrease, algorithm, forcing_function.
%   
%   nb                          Number of blocks.
%   maxfun                      Maximum of function evaluation.
%   maxfun_dim                  Factor of maximum of function evaluation regarding to dimenstions.
%   expand                      Expanding factor of step size.
%   shrink                      Shrinking factor of step size.
%   sufficient_decrease_factor  Factor of sufficient decrease condition.
%   StepTolerance               StepTolerance of step size. If step size is below StepTolerance,
%                               then the algorithm terminates.
%   ftarget                     If function value is below ftarget, then the algorithm terminates.
%   polling_inner               Polling strategy of each block.
%   with_memory                 Whether the cycling strategy employed in the opportunistic case
%                               memorizes the history or not.
%   cycling                     Cycling strategy employed in the opportunistic case.
%   accept_simple_decrease      Whether the algorithm accepts simple decrease to update xval and fval.
%   algorithm                   algorithm of BDS.
%   forcing_function            type of forcing function.
%
%   [XVAL, FVAL] = BDS(FUN, X0, OPTIONS) returns the value of the objective 
%   function FUN at the solution XVAL.
%
%   [XVAL, FVAL, EXITFLAG] = BDS(FUN, X0, OPTIONS) returns an EXITFLAG
%   that describes the exit condition. The information of EXITFLAG will be
%   given in output.message.
%
%   [XVAL, FVAL, EXITFLAG, OUTPUT] = BDS(FUN, X0, OPTIONS) returns a
%   structure OUTPUT with fields: fhist, xhist, alpha_hist, blocks_hist, funcCount.
%
%   fhist       History of function value
%   xhist       History of points that being computed
%   alpha_hist  History of step size every iteration
%   blocks_hist History of blocks that being visited
%   funcCount   The number of function evaluations.
%

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

% If fun is a string, then convert it to a function handle.
if ischarstr(fun)
    fun = str2func(fun);
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
end

% Get number of directions.
m = size(D, 2); 
% Get the number of blocks.
if isfield(options, "nb")
    nb = options.nb;
elseif strcmpi(options.Algorithm, "cbds") || strcmpi(options.Algorithm, "pbds")...
        || strcmpi(options.Algorithm, "rbds")
    % Default value is set as n, which is good for canonical with 2n directions. For
    % other situations, other value may be good.
    nb = n;
elseif strcmpi(options.Algorithm, "dspd") || strcmpi(options.Algorithm, "ds")
    nb = 1;
end

% Number of directions should be greater or equal to number of blocks.
nb = min(m, nb);
% Set indices of blocks as 1:nb.
block_indices = 1:nb;

% Set maxfun to the maximum number of function evaluations.
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

% Set value of expanding factor.
if isfield(options, "expand")
    expand = options.expand;
else
    expand = get_default_constant("expand");
end

% Set value of shrinking factor.
if isfield(options, "shrink")
    shrink = options.shrink;
else
    shrink = get_default_constant("shrink");
end

% Set value of sufficient decrease factor.
if isfield(options, "sufficient_decrease_factor")
    sufficient_decrease_factor = options.sufficient_decrease_factor;
else
    sufficient_decrease_factor = get_default_constant("sufficient_decrease_factor");
end

% Set type of forcing function.
if isfield(options, "forcing_function")
    forcing_function = options.forcing_function;
else
    forcing_function = get_default_constant("forcing_function");
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
% below StepTolerance, then the algorithm is stopped.
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

% Set the default value of polling_inner. This is the polling strategy
% employed within one block.
if ~isfield(options, "polling_inner")
    options.polling_inner = get_default_constant("polling_inner");
end

% Set the default value of cycling_inner. 
if isfield(options, "cycling_inner")
    cycling_inner = options.cycling_inner;
else
    cycling_inner = get_default_constant("cycling_inner");
end

% Set default value of shuffling_period. Default value of shuffling_period
% should be set as 1 since the algorithm visits all blocks for every iteration.
if strcmpi(options.Algorithm, "pbds") && isfield(options, "shuffling_period")
    shuffling_period = options.shuffling_period;
else
    shuffling_period = get_default_constant("shuffling_period");
end

% Set default value of replacement_delay. Default value of
% replacement_delay is set to be 0. Also, the value of replacement_delay
% should be less than or equal to nb-1, otherwise, there will not exist
% such block going to be visited after nb blocks have been visited.
if strcmpi(options.Algorithm, "rbds") && isfield(options, "replacement_delay")
    replacement_delay = min(options.replacement_delay, nb-1);
else
    replacement_delay = min(get_default_constant("replacement_delay"), nb-1);
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
% Initialize history of function values.
fhist = NaN(1, maxfun);
% Initialize history of points having been visited.
xhist = NaN(n, maxfun); 
% Initialize history of blocks having been visited.
block_hist = NaN(1, maxfun);
xval = x0; 
fval = eval_fun(fun, xval);
% Set number of function evaluations.
nf = 1; 
fhist(nf) = fval;
xhist(:, nf) = xval;

% Check whether ftarget is reached by fval. If this is the case, the
% computations afterwards need NOT be done.
if fval <= ftarget
    information = "FTARGET_REACHED";
    exitflag = get_exitflag(information);
    
    % The target function value has been reached at the very first function
    % evaluation. In this case, no further computation should be
    % entertained, and hence, no iteration should be run.
    maxit = 0;
end

% Start the actual computations.
% nb blocks have been explored after the number of iteration goes from k to k+1.
for iter = 1 : maxit
    % record the value of alpha_all of the current iteration in alpha_hist.
    alpha_hist(:, iter) = alpha_all;
    
    % Why iter-1? Since we will permute block_indices at the initial stage.
    if strcmpi(options.Algorithm, "pbds") && mod(iter - 1, shuffling_period) == 0
        % Make sure that `shuffling_period` is defined when `Algorithm` is "sbds".
        block_indices = randperm(nb);
    end
    
    % Get the block that are going to be visited in this iteration.
    if strcmpi(options.Algorithm, "rbds")
        % If replacement_delay is 0, then the algorithm will select a block randomly 
        % from block_indices for every iteration. If iter is equal to 1, then the block 
        % that we are going to visit is selected randomly from block_indices.
        if replacement_delay == 0 || iter == 1
            block_indices = randi([1, nb]);
        else
            % Record the number of blocks having been visited.
            num_visited = sum(~isnan(block_hist));
            % Get the number of blocks that we are going to exlude in the following selection.
            block_visited_slices_length = min(num_visited, replacement_delay);
            % Get the indices of blocks that we are going to exclude in the following selection.
            block_visited_slices = block_hist(num_visited-block_visited_slices_length+1:num_visited);
            % Set default value of initial block_indices.
            block_initial_indices = 1:nb;
            % Remove elements of block_indices appearing in block_visited_slice.
            block_real_indices = block_initial_indices(~ismember(block_initial_indices, block_visited_slices));
            % Generate a random index from block_real_indices.
            idx = randi(length(block_real_indices));
            block_indices = block_real_indices(idx);
        end
    end
    
    % Generate the searching set which its directions are uniformly distributed on the unit sphere for every iteration
    % when options.Algorithm is "dspd".
    if strcmpi(options.Algorithm, "dspd")
        rv = NaN(n, 1);
        for i = 1:n
             rv(i) = randn(1);
        end
        % Calculate the l2-norm of rv.
        norms = sqrt(sum(rv.^2, 1));
        % Normalize rv.
        rv = rv ./ norms;
        [Q, ~] = qr(rv);
        D = [Q, -Q];
    end

    for i = 1:length(block_indices)
        % In case of permutation.
        i_real = block_indices(i);
        
        % Recorder the number of blocks having been visited.
        num_visited = sum(~isnan(block_hist));
        % Update the block that going to be visited.
        block_hist(num_visited+1) = i_real;
        
        % Get indices in the i-th block.
        direction_indices = searching_set_indices{i_real}; 
        
        suboptions.maxfun = maxfun - nf;
        % Memory and cycling are needed since we may permutate indices in inner_direct_search.
        suboptions.cycling = cycling_inner;
        suboptions.with_memory = with_memory;
        suboptions.sufficient_decrease_factor = sufficient_decrease_factor;
        suboptions.ftarget = ftarget;
        suboptions.polling_inner = options.polling_inner;
        suboptions.accept_simple_decrease = accept_simple_decrease;
        suboptions.forcing_function = forcing_function;
        
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
    
end

% Set useful pieces on information about the solver's history in output.
output.funcCount = nf;
output.fhist = fhist(1:nf);
output.xhist = xhist(:, 1:nf);
output.alpha_hist = alpha_hist(:, 1:min(iter, maxit));
% Recorder the number of blocks visited.
num_blocks_visited = sum(~isnan(block_hist));
% Recorder the blocks visited.
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
