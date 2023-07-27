function [xval, fval, exitflag, output] = inner_direct_search(fun, ...
    xval, fval, D, direction_indices, alpha, options)
% inner_direct_search peforms a single iteration of traditional non-block
% direct search within a given block, for which the searching direction
% set is D.
%
% XVAL = INNER_DIRECT_SEARCH(FUN, XVAL, FVAL, D, ...
% DIRECTION_INDICES, ALPHA) attempts to find a XVAL to satifsy sufficient
% decrease condition.
%
% XVAL = INNER_DIRECT_SEARCH(FUN, XVAL, FVAL, D, ...
% DIRECTION_INDICES, ALPHA, OPTIONS) works with the structure OPTIONS.
% INNER_DIRECT_SEARCH uses these options: sufficient decrease factor,
% ftarget, polling, with_memory, cycling.
%
% [XVAL, FVAL] = INNER_DIRECT_SEARCH(...) returns the value of the
% objective function, described in FUN, at XVAL.
%
% [XVAL, FVAL, EXITFLAG] = INNER_DIRECT_SEARCH(...) returns an EXITFLAG
% that describes the exit condition. If it is normal, EXITFLAG will be NaN.
%
% [XVAL, FVAL, EXITFLAG, OUTPUT] = INNER_DIRECT_SEARCH(...) returns a
% structure OUTPUT with the number of function evaluations in OUTPUT.funcCount,
% the history of function evaluation in OUTPUT.fhist, the history of points
% in OUTPUT.xhist, boolean value of success, boolean value of terminate and
% direction_indices.
%
% success: success is initialized to be false. If success is updated to be true,
% it means that there at least exists some direction satisfying sufficient decrease.
%
% terminate: terminate is initialized to be false. If terminate is updated to be true,
% it means that either function evaluation is exhausted, or ftarget is reached.
%
% direction_indices: indices of directions of this block in D.


% Set options to an empty structure if it is not supplied.
if nargin < 7
    options = struct();
end

% Set the default sufficient decrease factor.
if isfield(options, "sufficient_decrease_factor")
    sufficient_decrease_factor = options.sufficient_decrease_factor;
else
    sufficient_decrease_factor = get_default_constant("sufficient_decrease_factor");
end

% Set the target on the objective function. If an evaluation of the
% objective function is below the target (the problem is unconstrained),
% then the algorithm is stopped.
if isfield(options, "ftarget")
    ftarget = options.ftarget;
else
    ftarget = get_default_constant("ftarget");
end

% Explain why NaN is good. It is possible that this function returns
% with exitflag=NaN and this is NOT a bug. This is because other situatons
% are corresponding to other normal values. Easy to see whether there is
% some bug related to exitflag.
exitflag = NaN;

% Initialize the computations.
n = length(xval);
num_directions = length(direction_indices);
fhist = NaN(1, num_directions);
xhist = NaN(n, num_directions);

% Initialize some parameters before entering the main computations.
nf = 0; % number of (inner) function evaluations
success = false; % the sufficient decrease condition is achieved
fbase = fval;
xbase = xval;
terminate = false;

% Start the main computations.
for j = 1 : num_directions
    % Stop the computations if no more function evaluation can be
    % performed. Note that this should be checked only before evaluating
    % the objective function.
    if nf >= options.maxfun
        terminate = true;
        exitflag = get_exitflag("MAXFUN_REACHED");
        break;
    end
    
    % Evaluate the objective function for the current polling direction.
    xnew = xbase+alpha*D(:, j);
    fnew = fun(xnew);
    nf = nf+1;
    fhist(nf) = fnew;
    xhist(:, nf) = xnew;
    
    % Stop the computations if the target value on the objective function
    % is achieved. Note that the comparison is done here because fnew may
    % be below ftarget without achieving a sufficient decrease.
    if fnew <= ftarget
        xval = xnew;
        fval = fnew;
        terminate = true;
        information = "FTARGET_REACHED";
        exitflag = get_exitflag(information);
        break;
    end
    
    % 1. Comment why the following line is wrong: if the following line is right, the complete
    %    polling will only receive success of the last direction, even if there exists success
    %    directions before it.
    % 2. What if we update fnew and xnew whenever there is a smple decrease?
    %success = (fnew <= fbase - sufficient_decrease_factor * alpha^2 / 2);
    
    sufficient_decrease = (fnew + sufficient_decrease_factor * alpha^2/2 < fbase);
    % Success is initialized to be false. Once there exists some direction satisfying sufficient
    % decrease, success will always be true.
    success = (success || sufficient_decrease);
    
    % For complete polling, fbase is fixed during all iterations in the block. It may happen that
    % where sufficient_decrease is true and fnew >= fval. For opportunistic polling, if
    % sufficient decrease is true, then the following points will not be explored.
    if (options.accept_simple_decrease || sufficient_decrease) && fnew < fval
        xval = xnew;
        fval = fnew;
    end
    
    % In the opportunistic case, if the current iteration achieves sufficient decrease,
    % stop the computations after cycling the indices of the polling
    % directions.
    if success && ~strcmpi(options.polling_inner, "complete")
        direction_indices = cycling(direction_indices, j, options.cycling, options.with_memory);
        break;
    end
end

% Set useful pieces on information about the history in output.
output.fhist = fhist(1:nf);
output.xhist = xhist(:, 1:nf);
output.nf = nf;
output.success = success;
output.direction_indices = direction_indices;
output.terminate = terminate;
end

function array = cycling(array, index, strategy, with_memory)
%CYCLING Permute an array according to different options
%   ARRAY = CYCLING(ARRAY, INDEX, STRATEGY, MEMORY) returns an array
%   that is a permutation of ARRAY according to INDEX, STRATEGY, and MEMORY.
%
%   ARRAY is the array to permute. It must be a vector.
%   INDEX is a number from -1, 1, ..., length(array). If INDEX = -1, then there is
%   no permutation.
%   MEMORY is a boolean value. If MEMORY is true, then the output ARRAY will
%   be obtained by permitting the ARRAY; otherwise, the input ARRAY will be
%   discarded and the output ARRAY will be obtained by permuting sort(ARRAY).
%   STRATEGY is a nonnegative integer from 0 to 4, indicating the strategy of the
%   permutation as follows.
%
%
%   0  No permutation.
%
%   1  The element of the index will be moved to the first element of array.
%
%   EXAMPLE
%   When array is 3 1 2 4 5, if index = 3, for with_memory situation,
%   array will be 2 3 1 4 5 after cycling; for nonwith_memory situaion, index
%   will be 2, sort(index) is 1 2 3 4 5 and array will be 2 1 3 4 5 after
%   cycling.
%
%   2  The element of the index and the following ones until end will be
%      moved ahead of array.
%
%   EXAMPLE
%   When array is 2 1 4 5 3, if index = 3, for with_memory situation,
%   array will be 4 5 3 2 1 after cycling; for nonwith_memory situaion, index
%   will be 4, sort(index) is 1 2 3 4 5 and array will be 4 5 1 2 3 after
%   cycling.
%
%   3  The element of the following ones after index until end will be
%      moved ahead of array.
%
%   EXAMPLE
%   When array is 2 1 4 5 3 and index = 3, for with_memory situation,
%   array will be 5 3 2 1 4 after cycling; for nonwith_memory situaion, index will
%   be 4, sort(index) is 1 2 3 4 5 and array will be 5 1 2 3 4 after cycling.
%
%   4  The element of the following one after index will be moved ahead of array.
%
%   EXAMPLE
%   array is 4 1 2 3 5, if index = 3, for with_memory situation, array will
%   be 3 4 1 2 5 after cycling; for nonwith_memory situaion, index will be 2,
%   sort(index) is 1 2 3 4 5 and array will be 3 1 2 4 5 after cycling.
%
%


% Precondition: If debug_flag is true, then pre-conditions is operated on
% input. If input_correctness is false, then assert may let the code crash.
debug_flag = is_debugging();
if debug_flag
    % Assert array is a real vector.
    [isrv, ~]  = isrealvector(array);
    assert(isrv);
    % Assert index is an integer.
    assert(isintegerscalar(index));
    % Assert strategy is a positive integer and less than or equal to 4.
    assert(isintegerscalar(strategy) && 0<=strategy && strategy<=4);
    % Assert with_memory is a boolean value.
    assert(islogicalscalar(with_memory));
end

%   If index < 0, then there is no "success_index" and there is no
%   permutation. If strategy == 0, then the permutation is unchanged.
if index < 0 || strategy == 0
    return;
end

% If with_memory is true, cycling_strategy will be operated on array. Otherwise,
% cycling_strategy will be operated on the array after sorting. In this case,
% the value of index will be the index corresponding to the array after sorting.
if ~with_memory
    [array, indices] = sort(array);
    index = find(indices == index);
end

switch strategy
    % If cycling_strategy is 1, the element of the index will be moved to
    % the first element of array. For example, if index = 3, array is 1 2 3 4 5,
    % then array will be 3 1 2 4 5 after cycling. For the case where
    % array is 3 1 2 4 5, if index = 3, for with_memory situation, array will
    % be 2 3 1 4 5 after cycling; for nonwith_memory situaion, index will
    % be 2 after executing the code paragraph above, sort(index)
    %is 1 2 3 4 5 and array will be 2 1 3 4 5 after cycling.
    case {1}
        array(1:index) = array([index, 1:index-1]);
        % If cycling_strategy is 2, the element of the index and the following
        % ones until end will be moved ahead of array. For example, if index = 3,
        % array is 1 2 3 4 5, then array will be 3 4 5 1 2 after cycling.
        % When array is 2 1 4 5 3, if index = 3, for with_memory
        % situation, array will be 4 5 3 2 1 after cycling; for nonwith_memory
        % situaion, index will be 4 after executing the paragraph above,
        % sort(index) is 1 2 3 4 5 and array will be 4 5 1 2 3 after cycling.
    case {2}
        array = array([index:end, 1:index-1]);
        % If cycling_strategy is 3, the element of the following ones after index
        % until end will be moved ahead of array. For example, if index = 3, array
        % is 1 2 3 4 5, then array will be 4 5 1 2 3 after cycling.
        % When array is 2 1 4 5 3 and index = 3, for with_memory
        % situation, array will be 5 3 2 1 4 after cycling; for nonwith_memory
        % situaion, index will be 4 after executing the paragraph above,
        % sort(index) is 1 2 3 4 5 and array will be 5 1 2 3 4 after cycling.
    case {3}
        array = array([index+1:end, 1:index]);
        % If cycling_strategy is 4, the element of the following one after index
        % will be moved ahead of array. For example, if index = 3, array
        % is 1 2 3 4 5, then array will be 4 1 2 3 5 after cycling.
        % For the case where array is 4 1 2 3 5, if index = 3, for with_memory
        % situation, array will be 3 4 1 2 5 after cycling; for nonwith_memory
        % situaion, index will be 2 after executing the paragraph above,
        % sort(index) is 1 2 3 4 5 and array will be 3 1 2 4 5 after cycling.
    case {4}
        if index ~= length(array)
            array(1:index+1) = array([index+1, 1:index]);
        end
end

if debug_flag
    % Assert array is a vector.
    [isrv, ~]  = isrealvector(array);
    assert(isrv);
end

end
