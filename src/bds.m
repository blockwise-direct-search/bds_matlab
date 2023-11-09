function [xopt, fopt, exitflag, output] = bds(fun, x0, options)
%BDS solves unconstrained optimization problems without using derivatives by blockwise direct search methods. 
%
%   BDS supports in MATLAB R2017b or later.
%   
%   XOPT = BDS(FUN, X0) returns an approximate minimizer XOPT of the function FUN, starting the 
%   calculations at X0. FUN must accept a vector input X and return a scalar.
%
%   XOPT = BDS(FUN, X0, OPTIONS) performs the computations with the options in OPTIONS. OPTIONS should be a
%   structure with the following fields.
%
%   Algorithm                   Algorithm to use. It can be "cbds" (cyclic blockwise direct search), 
%                               "pbds" (randomly permuted blockwise direct search), "rbds" (randomized 
%                               blockwise direct search), "ds" (the classical direct search without blocks).
%                               Default: "cbds".
%   nb                          Number of blocks. A positive integer. Default: n if Algorithm is "cbds", "pbds", 
%                               or "rbds", 1 if Algorithm is "ds".
%   maxfun                      Maximum of function evaluations. A positive integer. See also maxfun_factor.
%   maxfun_factor               Factor to define the maximum number of function evaluations as a multiple
%                               of the dimension of the problem. A positive integer. See also maxfun.
%                               The maximum of function evaluations is min(maxfun, maxfun_factor*n) if the user
%                               specify both maxfun and maxfun_factor; it is maxfun if the user only specifies 
%                               maxfun; it is maxfun_factor*n if the user only specifies maxfun_factor; it is
%                               min(get_default_constant("maxfun"), get_default_constant("maxfun_factor")*n) if 
%                               the user specifies neither maxfun nor maxfun_factor.
%   direction_set               A set of directions used for polling. It can be "canonical", "identity", or
%                               a matrix of n rows. See get_direction_set.m for details. Default: "canonical".
%   expand                      Expanding factor of step size. A real number no less than 1. Default: 2.
%   shrink                      Shrinking factor of step size. A positive number less than 1. Default: 0.5.
%   forcing_function            The forcing function used for deciding whether the step achieves a sufficient
%                               decrease. A function handle. Default: @(alpha) alpha^2. See also reduction_factor. 
%   reduction_factor            Factors multiplied to the forcing function when deciding whether the step achieves
%                               a sufficient decrease. A 3-dimentional vector such that 
%                               reduction_factor(1) <= reduction_factor(2) <= reduction_factor(3),
%                               reduction_factor(1) >= 0, and reduction_factor(2) > 0.
%                               reduction_factor(0) is used for deciding whether to update the base point; 
%                               reduction_factor(1) is used for deciding whether to shrink the step size; 
%                               reduction_factor(2) is used for deciding whether to expand the step size.
%                               Default: [0, eps, eps]. See also forcing_function.
%   StepTolerance               Lower bound of the step size. If the step size is smaller than StepTolerance,
%                               then the algorithm terminates. A (small) positive number. Default: 1e-10.
%   ftarget                     Target of the function value. If the function value is smaller than or equal to
%                               ftarget, then the algorithm terminates. A real number. Default: -Inf.
%   polling_inner               Polling strategy in each block. It can be "complete" or "opportunistic". 
%                               Default: "opportunistic".
%   cycling_inner               Cycling strategy employed within each block. It is used only when polling_inner 
%                               is "opportunistic". It can be 0, 1, 2, 3, 4. See cycling.m for details. 
%                               Default: 3.
%   with_cycling_memory         Whether the cycling strategy within each block memorizes the history or not. 
%                               It is used only when polling_inner is "opportunistic". Default: true.
%   shuffling_period            It is only used in PBDS, which shuffles the blocks every shuffling_period 
%                               iterations. A positive integer. Default: 1.   
%   replacement_delay           It is only used for RBDS. Suppose that replacement_delay is r. If block i
%                               is selected at iteration k, then it will not be selected at iterations 
%                               k+1, ..., k+r. An integer between 0 and nb-1. Default: 0.
%   seed                        The seed for permuting blocks in PBDS or randomly choosing one block in RBDS.
%                               It is only for reproducibility in experiments. A positive integer.
%   output_xhist                Whether the history of points visited is returned or not. Default: false.
%   output_alpha_hist           Whether the history of step sizes is returned or not. Default: false.
%   output_block_hist           Whether the history of blocks visited is returned or not. Default: false.
%
%   [XOPT, FOPT] = BDS(...) also returns the value of the objective function FUN at the 
%   solution XOPT.
%
%   [XOPT, FOPT, EXITFLAG] = BDS(...) returns an EXITFLAG that describes the exit 
%   condition. The possible values of EXITFLAG are 0, 1, 2, 3.
%
%   0    The StepTolerance of the step size is reached.
%   1    The maximum number of function evaluations is reached.
%   2    The target of the objective function is reached.
%   3    The maximum number of iterations is reached.
%
%   [XOPT, FOPT, EXITFLAG, OUTPUT] = BDS(...) returns a
%   structure OUTPUT with the following fields: 
%
%   fhist        History of function values.
%   xhist        History of points visited.
%   alpha_hist   History of step size for every iteration.
%   blocks_hist  History of blocks visited.
%   funcCount    The number of function evaluations.
%   message      The information of EXITFLAG.
%
%   Copyright 2023, Haitian Li and Zaikun Zhang.
%   All rights reserved.
%

% Set options to an empty structure if it is not provided.
if nargin < 3
    options = struct();
end

% Transpose x0 if it is a row.
x0_is_row = isrow(x0);
x0 = double(x0(:));

% Check the inputs of the user when debug_flag is true.
debug_flag = is_debugging();
if debug_flag
    verify_preconditions(fun, x0, options);
end

% If FUN is a string, then convert it to a function handle.
if ischarstr(fun)
    fun = str2func(fun);
end
% Redefine fun to accept columns if x0 is a row, as we use columns internally.
if x0_is_row
    fun = @(x)fun(x');
end

% Set the polling directions in D.
n = length(x0);

if ~isfield(options, "Algorithm")
    options.Algorithm = get_default_constant("Algorithm");
end

% Get the direction set of directions.
D = get_direction_set(n, options);

% Set the value of expanding factor.
if isfield(options, "expand")
    expand = options.expand;
else
    expand = get_default_constant("expand");
end

% Set the value of shrinking factor.
if isfield(options, "shrink")
    shrink = options.shrink;
else
    shrink = get_default_constant("shrink");
end

% Get the number of directions.
m = size(D, 2);
 
% Get the number of blocks.
if isfield(options, "nb")
    % The number of directions should be greater or equal to the number of blocks.
    nb = min(m, options.nb);
elseif strcmpi(options.Algorithm, "cbds") || strcmpi(options.Algorithm, "pbds")...
        || strcmpi(options.Algorithm, "rbds")
    % Default value is set as n, which is good for canonical with 2n directions. For
    % other situations, other value may be good.
    nb = n;
elseif strcmpi(options.Algorithm, "ds")
    nb = 1;
end

% Set indices of blocks as 1:nb.
block_indices = 1:nb;

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

% Set the value of reduction factor.
if isfield(options, "reduction_factor")
    reduction_factor = options.reduction_factor;
else
    reduction_factor = get_default_constant("reduction_factor");
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

% Set the value of shuffling_period, which shuffles the blocks every shuffling_period iterations.
if strcmpi(options.Algorithm, "pbds") && isfield(options, "shuffling_period")
    shuffling_period = options.shuffling_period;
else
    shuffling_period = get_default_constant("shuffle_period");
end

% Set the value of replacement_delay. The default value of replacement_delay is set to 0. 
if strcmpi(options.Algorithm, "rbds") && isfield(options, "replacement_delay")
    replacement_delay = min(options.replacement_delay, nb-1);
else
    replacement_delay = min(get_default_constant("replacement_delay"), nb-1);
end

% Set the boolean value of WITH_CYCLING_MEMORY. 
% WITH_CYCLING_MEMORY is only used when we need to permute the directions_indices. If
% WITH_CYCLING_MEMORY is true, then we will permute the directions_indices by using the
% directions_indices of the previous iteration. Otherwise, we will permute the directions_indices
% with the initial directions_indices of ascending orders.
if isfield(options, "with_cycling_memory")
    with_cycling_memory = options.with_cycling_memory;
else
    with_cycling_memory = get_default_constant("with_cycling_memory");
end

% Initialize the step sizes and alpha_hist, which is the history of step sizes.
if isfield(options, "output_alpha_hist")
    output_alpha_hist = options.output_alpha_hist;
else
    output_alpha_hist = get_default_constant("output_alpha_hist");
end

if output_alpha_hist
    try
        alpha_hist = NaN(nb, maxit);
    catch
        output_alpha_hist = false;
        warning("The size of alpha_hist exceeds the maximum of memory size limit.")
    end
end

if isfield(options, "alpha_init")
    if length(options.alpha_init) == 1
        alpha_all = options.alpha_init*ones(nb, 1);
    elseif length(options.alpha_init) == nb
        alpha_all = options.alpha_init;
    else
        error("The length of alpha_init should be equal to nb or equal to 1.");
    end
    % Try alpha_all = 0.5 * max(abs(x0), 1) in the canonical case.
elseif isfield(options, "alpha_init_scaling") && options.alpha_init_scaling
    %alpha_all = 0.1 * ones(nb, 1);
    %alpha_all(x0 ~= 0) = 0.1 * abs(x0(x0 ~= 0));
    % alpha_all = 0.5 * max(abs(x0), ones(nb, 1));
    alpha_all = 0.1 * max(1e-3, abs(x0));
else
    alpha_all = ones(nb, 1);
end

% Decide which polling direction belongs to which block.
direction_set_indices = divide_direction_set(m, nb);

% Initialize the history of function values.
fhist = NaN(1, maxfun);

% Initialize the history of points visited.
if isfield(options, "output_xhist")
    output_xhist = options.output_xhist;
else
    output_xhist = get_default_constant("output_xhist");
end

if output_xhist
    try
        xhist = NaN(n, maxfun);
    catch
        output_xhist = false;
        warning("xhist will be not included in the output due to the limit of memory.");
    end
end

if isfield(options, "output_block_hist")
    output_block_hist = options.output_block_hist;
else
    output_block_hist = get_default_constant("output_block_hist");
end

% Initialize the history of blocks visited. If output_block_hist is true, then we will return the
% history of blocks visited in the output.
block_hist = NaN(1, maxfun);

% To avoid that the users bring some randomized strings.
if ~isfield(options, "seed")
    options.seed = get_default_constant("seed");
end
random_stream = RandStream("mt19937ar", "Seed", options.seed);

% Initialize the exitflag where the maximum number of iterations is reached. 
exitflag = get_exitflag("MAXIT_REACHED");
xbase = x0; 
fbase = eval_fun(fun, xbase);
% Set the number of function evaluations.
nf = 1; 
if output_xhist
    xhist(:, nf) = xbase;
end
fhist(nf) = fbase;
xopt = xbase;
fopt = fbase;
terminate = false;

% Check whether FTARGET is reached by fopt. If it is true, then terminate.
if fopt <= ftarget
    information = "FTARGET_REACHED";
    exitflag = get_exitflag(information);
    
    % FTARGET has been reached at the very first function evaluation. 
    % In this case, no further computation should be entertained, and hence, 
    % no iteration should be run.
    maxit = 0;
end

for iter = 1:maxit
    % Record the value of alpha_all of the current iteration in alpha_hist.
    if output_alpha_hist
        alpha_hist(:, iter) = alpha_all;
    end
    
    % Shuffle the blocks every shuffling_period iterations.
    % Why iter-1? Since we will permute block_indices at the initial stage.
    if strcmpi(options.Algorithm, "pbds") && mod(iter - 1, shuffling_period) == 0
        % Make sure that shuffling_period is defined when the Algorithm is "pbds".
        block_indices = random_stream.randperm(nb);
    end
    
    % Get the block that is going to be visited.
    if strcmpi(options.Algorithm, "rbds")
        % If replacement_delay is 0, then select a block randomly from block_indices for 
        % each iteration. If iter is equal to 1, then the block that we are going to visit
        % is selected randomly from block_indices.
        if replacement_delay == 0 || iter == 1
            block_indices = random_stream.randi([1, nb]);
        else
            % Record the number of blocks visited.
            num_visited = sum(~isnan(block_hist));
            % Get the number of blocks that we are going to exclude in the following selection.
            block_visited_slices_length = min(num_visited, replacement_delay);
            % Get the indices of blocks that we are going to exclude in the following selection.
            block_visited_slices = block_hist(num_visited-block_visited_slices_length+1:num_visited);
            % Set the default value of initial block_indices.
            block_initial_indices = 1:nb;
            % Remove elements of block_indices appearing in block_visited_slice.
            block_real_indices = block_initial_indices(~ismember(block_initial_indices, block_visited_slices));
            % Generate a random index from block_real_indices.
            idx = random_stream.randi(length(block_real_indices));
            block_indices = block_real_indices(idx);
        end
    end
    
    for i = 1:length(block_indices)
        % If block_indices is 1 3 2, then block_indices(2) = 3, which is the real block that we are
        % going to visit.
        i_real = block_indices(i);
        
        % Record the number of blocks visited.
        num_visited = sum(~isnan(block_hist));
        % Record the block that is going to be visited.
        block_hist(num_visited+1) = i_real;       
       
        % Get indices of directions in the i-th block.
        direction_indices = direction_set_indices{i_real}; 
        
        suboptions.maxfun = maxfun - nf;
        suboptions.cycling_inner = cycling_inner;
        suboptions.with_cycling_memory = with_cycling_memory;
        suboptions.reduction_factor = reduction_factor;
        suboptions.forcing_function = forcing_function;
        suboptions.ftarget = ftarget;
        suboptions.polling_inner = options.polling_inner;
        
        [sub_xopt, sub_fopt, sub_exitflag, sub_output] = inner_direct_search(fun, xbase,...
            fbase, D(:, direction_indices), direction_indices,...
            alpha_all(i_real), suboptions);


        % Update the history of step size.
        if output_alpha_hist
            alpha_hist(:, iter) = alpha_all;
        end
        
        % Store the history of the evaluations by inner_direct_search, 
        % and accumulate the number of function evaluations.
        fhist((nf+1):(nf+sub_output.nf)) = sub_output.fhist;
        if output_xhist
            xhist(:, (nf+1):(nf+sub_output.nf)) = sub_output.xhist;
        end

        nf = nf+sub_output.nf;
        
        % Update the step sizes and store the history of step sizes.
        % fbase and xbase are used for the computation in the block. fopt and
        % xopt are always the best function value and point so far.
        if sub_fopt + reduction_factor(3) * forcing_function(alpha_all(i_real)) < fbase
            alpha_all(i_real) = expand * alpha_all(i_real);
        elseif sub_fopt + reduction_factor(2) * forcing_function(alpha_all(i_real)) >= fbase
            alpha_all(i_real) = shrink * alpha_all(i_real);
        end

        if (reduction_factor(1) <= 0 && sub_fopt < fbase) || sub_fopt + reduction_factor(1) * forcing_function(alpha_all(i_real)) < fbase
            xbase = sub_xopt;
            fbase = sub_fopt;
        end

        % Update xopt and fopt.
        if sub_fopt < fopt
            xopt = sub_xopt;
            fopt = sub_fopt;
        end
                        
        % Retrieve the order of the polling directions in the direction set.
        direction_set_indices{i_real} = sub_output.direction_indices;

        % If sub_output.terminate is true, then inner_direct_search returns 
        % boolean value of terminate because either the maximum number of function
        % evaluations or the target of the objective function value is reached. 
        % In both cases, the exitflag is set by inner_direct_search.
        if sub_output.terminate
            terminate = true;
            exitflag = sub_exitflag;
            break;
        end

        % Terminate the computations if the largest component of step size is below a
        % given StepTolerance.
        if max(alpha_all) < alpha_tol
            terminate = true;
            exitflag = get_exitflag("SMALL_ALPHA");
            break;
        end 
    end
    
    % Check whether one of SMALL_ALPHA, MAXFUN_REACHED, and FTARGET_REACHED is reached.
    if terminate
        break;
    end
    
end

% Truncate HISTORY into a vector of nf length.
output.funcCount = nf;
output.fhist = fhist(1:nf);

if output_xhist
    output.xhist = xhist(:, 1:nf);
end

if output_alpha_hist
    output.alpha_hist = alpha_hist(1:min(iter, maxit));
end

num_blocks_visited = sum(~isnan(block_hist));
if output_block_hist
    output.blocks_hist = block_hist(1:num_blocks_visited);
end

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

% Transpose xopt if x0 is a row.
if x0_is_row
    xopt = xopt';
end

% verify_postconditions is to detect whether the output is in the right form when debug_flag is true.
if debug_flag
    verify_postconditions(fun, xopt, fopt, exitflag, output);
end
