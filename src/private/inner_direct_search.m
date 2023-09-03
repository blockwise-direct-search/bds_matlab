function [xval, fval, exitflag, output] = inner_direct_search(fun, ...
    xval, fval, D, direction_indices, alpha, options)
%INNER_DIRECT_SEARCH peforms a single iteration of classical direct search 
%   within a given block.
%
%   XVAL = INNER_DIRECT_SEARCH(FUN, XVAL, FVAL, D, DIRECTION_INDICES, ALPHA, OPTIONS)
%   returns a XVAL.
%
%   XVAL = INNER_DIRECT_SEARCH(FUN, XVAL, FVAL, D, ...
%   DIRECTION_INDICES, ALPHA, OPTIONS) works with the structure OPTIONS, which includes
%   sufficient_decrease_factor, ftarget, polling, with_memory, cycling, forcing_function.
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
%   TERMINATE is initialized to be false. If terminate is updated to be true,
%   it means that either the number of function evaluations reaches maxfun, or ftarget is reached.
%
%   DIRECTION_INDICES is indices of directions of this block in D.


% Set options to an empty structure if it is not provided.
if nargin < 7
    options = struct();
end

% Set the value of sufficient decrease factor.
if isfield(options, "sufficient_decrease_factor")
    sufficient_decrease_factor = options.sufficient_decrease_factor;
else
    sufficient_decrease_factor = get_default_constant("sufficient_decrease_factor");
end

% Set the type of forcing function.
if isfield(options, "forcing_function")
    forcing_function = options.forcing_function;
else
    forcing_function = get_default_constant("forcing_function");
end

% Set ftarget of objective function.
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

% Initialize some parameters before entering the loop.
n = length(xval);
num_directions = length(direction_indices);
fhist = NaN(1, num_directions);
xhist = NaN(n, num_directions);
nf = 0; 
success = false; 
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
    if strcmpi(forcing_function, "quadratic")
        % In this case, forcing function is c*alpha^2/2.
        sufficient_decrease = (fnew + sufficient_decrease_factor * alpha^2/2 < fbase);
        % In this case, forcing function is 0.
    elseif strcmpi(forcing_function, "zero")
        sufficient_decrease = (fnew < fbase);
    end    

    % Success is initialized to be false. Once there exists some direction satisfying sufficient
    % decrease, success will always be true inside this function.
    success = (success || sufficient_decrease);
    
    % If options.accept_simple_decrease is true, then we will accept xnew
    % and fnew as xval and fval respectively as long as fnew < fval. Otherwise,
    % we will only accept xnew and fnew when they meet both sufficient decrease and 
    % fnew < fval simultaneously. For complete polling, fbase is fixed during all 
    % iterations in the block. So there are some cases where sufficient_decrease is true 
    % and fnew >= fval. For opportunistic polling, as long as sufficient decrease is true, 
    % then the remaining polling points will not be explored. 
    if (options.accept_simple_decrease || sufficient_decrease) && fnew < fval
        xval = xnew;
        fval = fnew;
    end
    
    % In the opportunistic case, if the current iteration achieves sufficient decrease,
    % stop the computations after cycling the indices of the polling directions.
    if sufficient_decrease && ~strcmpi(options.polling_inner, "complete")
        direction_indices = cycling(direction_indices, j, options.cycling, options.with_memory);
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
end

function array = cycling(array, index, strategy, with_memory)
%CYCLING permutes an array according to different options.
%   ARRAY = CYCLING(ARRAY, INDEX, STRATEGY, MEMORY) returns an array
%   that is a permutation of ARRAY according to INDEX, STRATEGY, and MEMORY.
%
%   ARRAY is the array to permute. It must be a vector.
%   INDEX is a number from -1, 1, 2, ..., length(array). If INDEX = -1, then there is
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

% Check whether input is given in correct type when debug_flag is true. 
debug_flag = is_debugging();
if debug_flag
    % Array should be a real vector.
    [isrv, ~]  = isrealvector(array);
    if ~isrv
        error("Array is not a real vector.");
    end
    % Index should be an integer.
    if ~isintegerscalar(index)
        error("Index is not an integer.");
    end
    % Strategy should be a positive integer and less than or equal to 4.
    if ~isintegerscalar(strategy) || strategy < 0 || strategy > 4
        error("Strategy is not a positive integer or less than or equal to 4.");
    end
    % With_memory should be boolean value.
    if ~islogicalscalar(with_memory)
        error("With_memory is not a boolean value.");
    end
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

% Check whether ARRAY is a vector or not when debug_flag is true.
if debug_flag
    % Array should be a vector.
    [isrv, ~]  = isrealvector(array);
    if ~isrv
        error("Array is not a real vector.");
    end
end

end
