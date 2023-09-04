function array = cycling(array, index, strategy, with_cycling_memory)
%   CYCLING Permute an array according to different options
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
%   0  No permutation. 
%
%   1  The element of the index will be moved to the first element of array.
%
%   EXAMPLE 
%   When array is 3 1 2 4 5, if index = 3, for with_cycling_memory situation, 
%   array will be 2 3 1 4 5 after cycling; for nonwith_cycling_memory situaion, index 
%   will be 2, sort(index) is 1 2 3 4 5 and array will be 2 1 3 4 5 after 
%   cycling.
%
%   2  The element of the index and the following ones until end will be
%      moved ahead of array.
%
%   EXAMPLE
%   When array is 2 1 4 5 3, if index = 3, for with_cycling_memory situation, 
%   array will be 4 5 3 2 1 after cycling; for nonwith_cycling_memory situaion, index 
%   will be 4, sort(index) is 1 2 3 4 5 and array will be 4 5 1 2 3 after 
%   cycling.
%
%   3  The element of the following ones after index until end will be 
%      moved ahead of array.
%   
%   EXAMPLE 
%   When array is 2 1 4 5 3 and index = 3, for with_cycling_memory situation, 
%   array will be 5 3 2 1 4 after cycling; for nonwith_cycling_memory situaion, index will
%   be 4, sort(index) is 1 2 3 4 5 and array will be 5 1 2 3 4 after cycling.
%   
%   4  The element of the following one after index will be moved ahead of array.
%   
%   EXAMPLE  
%   array is 4 1 2 3 5, if index = 3, for with_cycling_memory situation, array will 
%   be 3 4 1 2 5 after cycling; for nonwith_cycling_memory situaion, index will be 2, 
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
    % Assert with_cycling_memory is a boolean value.
    assert(islogicalscalar(with_cycling_memory));
end

%   If index < 0, then there is no "success_index" and there is no
%   permutation. If strategy == 0, then the permutation is unchanged.
if index < 0 || strategy == 0
    return;
end

% If with_cycling_memory is true, cycling_strategy will be operated on array. Otherwise,
% cycling_strategy will be operated on the array after sorting. In this case,
% the value of index will be the index corresponding to the array after sorting.
if ~with_cycling_memory
    [array, indices] = sort(array);
    index = find(indices == index);
end

switch strategy
    % If cycling_strategy is 1, the element of the index will be moved to 
    % the first element of array. For example, if index = 3, array is 1 2 3 4 5, 
    % then array will be 3 1 2 4 5 after cycling. For the case where 
    % array is 3 1 2 4 5, if index = 3, for with_cycling_memory situation, array will 
    % be 2 3 1 4 5 after cycling; for nonwith_cycling_memory situaion, index will 
    % be 2 after executing the code paragraph above, sort(index) 
    %is 1 2 3 4 5 and array will be 2 1 3 4 5 after cycling.
    case {1}
        array(1:index) = array([index, 1:index-1]);
    % If cycling_strategy is 2, the element of the index and the following 
    % ones until end will be moved ahead of array. For example, if index = 3, 
    % array is 1 2 3 4 5, then array will be 3 4 5 1 2 after cycling. 
    % When array is 2 1 4 5 3, if index = 3, for with_cycling_memory 
    % situation, array will be 4 5 3 2 1 after cycling; for nonwith_cycling_memory 
    % situaion, index will be 4 after executing the paragraph above, 
    % sort(index) is 1 2 3 4 5 and array will be 4 5 1 2 3 after cycling.
    case {2}
        array = array([index:end, 1:index-1]);
    % If cycling_strategy is 3, the element of the following ones after index
    % until end will be moved ahead of array. For example, if index = 3, array 
    % is 1 2 3 4 5, then array will be 4 5 1 2 3 after cycling. 
    % When array is 2 1 4 5 3 and index = 3, for with_cycling_memory 
    % situation, array will be 5 3 2 1 4 after cycling; for nonwith_cycling_memory 
    % situaion, index will be 4 after executing the paragraph above, 
    % sort(index) is 1 2 3 4 5 and array will be 5 1 2 3 4 after cycling.
    case {3}
        array = array([index+1:end, 1:index]);
    % If cycling_strategy is 4, the element of the following one after index
    % will be moved ahead of array. For example, if index = 3, array 
    % is 1 2 3 4 5, then array will be 4 1 2 3 5 after cycling. 
    % For the case where array is 4 1 2 3 5, if index = 3, for with_cycling_memory 
    % situation, array will be 3 4 1 2 5 after cycling; for nonwith_cycling_memory 
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