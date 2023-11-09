function array = cycling(array, index, strategy, with_cycling_memory)
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
    %   1  The element of the index will be moved to the first element of the array.
    %
    %   EXAMPLE
    %   When the array is a3 a1 a2 a4 a5, if index = 3, the array will be a2 a3 a1 a4 a5 after cycling 
    %   when with_cycling_memory is true; if index = 2, sort(index) is a1 a2 a3 a4 a5 
    %   and array will be a2 a1 a3 a4 a5 after cycling when with_cycling_memory is false.
    %
    %   2  The element of the index and the following ones until the end will be
    %      moved ahead of the array.
    %
    %   EXAMPLE
    %   When the array is a2 a1 a4 a5 a3, if index = 3, the array will be a4 a5 a3 a2 a1 after cycling 
    %   when with_cycling_memory is true; if index = 4, sort(index) is a1 a2 a3 a4 a5 
    %   and array will be a4 a5 a1 a2 a3 after cycling when with_cycling_memory is false. 
    %
    %   3  The element of the following ones after the index until the end will be
    %      moved ahead of the array.
    %
    %   EXAMPLE
    %   When the array is a2 a1 a4 a5 a3 and index = 3, the array will be a5 a3 a2 a1 a4 after cycling 
    %   when with_cycling_memory is true; if index = 4, sort(index) is a1 a2 a3 a4 a5 
    %   and array will be a5 a1 a2 a3 a4 after cycling when with_cycling_memory is false.
    %
    %   4  The element of the following one after the index will be moved ahead of the array.
    %
    %   EXAMPLE
    %   array is a4 a1 a2 a3 a5, if index = 3, the array will be a3 a4 a1 a2 a5 after cycling when 
    %   with_cycling_memory is true; if index = 2, sort(index) is a1 a2 a3 a4 a5 and 
    %   array will be a3 a1 a2 a4 a5 after cycling when with_cycling_memory is false.
    %
    
    % Check whether the input is given in the correct type when debug_flag is true. 
    debug_flag = is_debugging();
    if debug_flag
        % Array should be a real vector.
        if ~isrealvector(array)
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
        % With_memory should be a boolean value.
        if ~islogicalscalar(with_cycling_memory)
            error("With_memory is not a boolean value.");
        end
    end
    
    %   If index < 0, then there is no "success_index" and there is no
    %   permutation. If strategy == 0, then the permutation is unchanged.
    if index < 0 || strategy == 0
        return;
    end
    
    % If with_cycling_memory is true, cycling_strategy will be operated on array. Otherwise,
    % cycling_strategy will be operated on the array after sorting. In this case,
    % the value of the index will be the index corresponding to the array after sorting.
    if ~with_cycling_memory
        [array, indices] = sort(array);
        index = find(indices == index);
    end
    
    switch strategy
        % If cycling_strategy is 1, the element of the index will be moved to
        % the first element of the array. For example, if index = 3, array is a1 a2 a3 a4 a5,
        % then the array will be a3 a1 a2 a4 a5 after cycling. For the case where
        % array is a3 a1 a2 a4 a5, if index = 3, the array will be a2 a3 a1 a4 a5 after cycling 
        % when with_cycling_memory is true; when with_cycling_memory is false, index will
        % be 2 after executing the code paragraph above, sort(index)
        %is a1 a2 a3 a4 a5 and the array will be a2 a1 a3 a4 a5 after cycling.
        case {1}
            array(1:index) = array([index, 1:index-1]);
            % If cycling_strategy is 2, the element of the index and the following
            % ones until the end will be moved ahead of the array. For example, if index = 3,
            % array is a1 a2 a3 a4 a5, then the array will be a3 a4 a5 a1 a2 after cycling.
            % When the array is a2 a1 a4 a5 a3, if index = 3, the array will be a4 a5 a3 a2 a1 after cycling
            % when with_cycling_memory is true; when with_cycling_memory
            % is false, the index will be 4 after executing the paragraph above,
            % sort(index) is a1 a2 a3 a4 a5 and the array will be a4 a5 a1 a2 a3 after cycling.
        case {2}
            array = array([index:end, 1:index-1]);
            % If cycling_strategy is 3, the element of the following ones after the index
            % until the end will be moved ahead of the array. For example, if index = 3, array
            % is a1 a2 a3 a4 a5, then the array will be a4 a5 a1 a2 a3 after cycling.
            % When the array is a2 a1 a4 a5 a3 and index = 3, the array will be a5 a3 a2 a1 a4 after cycling
            % when with_cycling_memory is true; when with_cycling_memory is false,
            % index will be 4 after executing the paragraph above,
            % sort(index) is a1 a2 a3 a4 a5 and the array will be a5 a1 a2 a3 a4 after cycling.
        case {3}
            array = array([index+1:end, 1:index]);
            % If cycling_strategy is 4, the element of the following one after the index
            % will be moved ahead of the array. For example, if index = 3, array
            % is a1 a2 a3 a4 a5, then the array will be a4 a1 a2 a3 a5 after cycling.
            % For the case where the array is a4 a1 a2 a3 a5, if index = 3, the array will 
            % be a3 a4 a1 a2 a5 after cycling when with_cycling_memory is true; when 
            % with_cycling_memory is false, index will be 2 after executing the paragraph above,
            % sort(index) is a1 a2 a3 a4 a5 and the array will be a3 a1 a2 a4 a5 after cycling.
        case {4}
            if index ~= length(array)
                array(1:index+1) = array([index+1, 1:index]);
            end
    
    % Check whether ARRAY is a vector or not when debug_flag is true.
    if debug_flag
        % Array should be a vector.
        if ~isrealvector(array)
            error("Array is not a real vector.");
        end
    end
    
    end