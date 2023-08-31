function tests = bds_unit_test
    %BDS_TEST runs all the test functions in this file.
    %   To run these tests, simply type "runtests" in the command window. To
    %   create a new test function in this file with a name that starts or
    %   finishes with "test" (case insensitive). For more info, see
    %
    %       https://www.mathworks.com/help/matlab/matlab_prog/write-simple-test-case-with-functions.html
    
    fullpath = mfilename("fullpath"); 
    [path, ~] = fileparts(fullpath);
    cd(path);
    cd ../src/private
    
    tests = functiontests(localfunctions);
    
    end
    
    function cycling_test(testCase)
    %CYCLING_TEST tests the file private/cycling.m
    
    % The following must not cycle the array.
    array = [1, 2, 3, 4, 5];
    for memory = [true, false]
        for strategy = 0:4
            verifyEqual(testCase, cycling(array, -1, strategy, memory), array)
        end
        for index = 1:length(array)
            verifyEqual(testCase, cycling(array, index, 0, memory), array)
        end
    end
    
    % The following are the ones commented in cycling.m
    array = [1, 2, 3, 4, 5];
    for memory = [true, false]
        verifyEqual(testCase, cycling(array, 3, 1, memory), [3, 1, 2, 4, 5])
        verifyEqual(testCase, cycling(array, 3, 2, memory), [3, 4, 5, 1, 2])
        verifyEqual(testCase, cycling(array, 3, 3, memory), [4, 5, 1, 2, 3])
        verifyEqual(testCase, cycling(array, 3, 4, memory), [4, 1, 2, 3, 5])
    end
    
    % The following tests the memory argument
    % Interestingly, these tests show that the comment in cycling.m stating
    % that the strategy 2 and 3 are not influenced by the value of memory is
    % wrong!!
    array = [2, 1, 4, 5, 3];
    verifyEqual(testCase, cycling(array, 3, 1, true), [4, 2, 1, 5, 3])
    verifyEqual(testCase, cycling(array, 3, 1, false), [4, 1, 2, 3, 5])
    verifyEqual(testCase, cycling(array, 3, 2, true), [4, 5, 3, 2, 1])
    verifyEqual(testCase, cycling(array, 3, 2, false), [4, 5, 1, 2, 3])
    verifyEqual(testCase, cycling(array, 3, 3, true), [5, 3, 2, 1, 4])
    verifyEqual(testCase, cycling(array, 3, 3, false), [5, 1, 2, 3, 4])
    
    % The following tests the strategy 5
    array = [1, 2, 3, 4, 5, 4, 3, 2, 1];
    
    % The following are the ones commented in cycling.m
    for memory = [true, false]
        verifyEqual(testCase, cycling(array, 1, 5, memory), [2, 3, 4, 5, 4, 3, 2, 1, 1]);
        verifyEqual(testCase, cycling(array, 2, 5, memory), [3, 4, 5, 4, 3, 2, 1, 1, 2]);
        verifyEqual(testCase, cycling(array, 3, 5, memory), [4, 5, 4, 3, 2, 1, 1, 2, 3]);
        verifyEqual(testCase, cycling(array, 4, 5, memory), [5, 4, 3, 2, 1, 1, 2, 3, 4]);
        verifyEqual(testCase, cycling(array, 5, 5, memory), [4, 3, 2, 1, 1, 2, 3, 4, 5]);
        verifyEqual(testCase, cycling(array, 6, 5, memory), [3, 2, 1, 1, 2, 3, 4, 5, 4]);
        verifyEqual(testCase, cycling(array, 7, 5, memory), [2, 1, 1, 2, 3, 4, 5, 4, 3]);
        verifyEqual(testCase, cycling(array, 8, 5, memory), [1, 1, 2, 3, 4, 5, 4, 3, 2]);
        verifyEqual(testCase, cycling(array, 9, 5, memory), [1, 2, 3, 4, 5, 4, 3, 2, 1]);
    end
    
    array = [3, 2, 2, 1, 1];
    verifyEqual(testCase, cycling(array, 1, 5, false), [2, 1, 1, 2, 3]);
    verifyEqual(testCase, cycling(array, 2, 5, false), [3, 2, 1, 1, 2]);
    verifyEqual(testCase, cycling(array, 3, 5, false), [1, 1, 2, 3, 2]);
    verifyEqual(testCase, cycling(array, 4, 5, false), [2, 3, 2, 1, 1]);
    verifyEqual(testCase, cycling(array, 5, 5, false), [1, 2, 3, 2, 1]);
    verifyEqual(testCase, cycling(array, 1, 5, true), [2, 2, 1, 1, 3]);
    verifyEqual(testCase, cycling(array, 2, 5, true), [2, 1, 1, 3, 2]);
    verifyEqual(testCase, cycling(array, 3, 5, true), [1, 1, 3, 2, 2]);
    verifyEqual(testCase, cycling(array, 4, 5, true), [1, 3, 2, 2, 1]);
    verifyEqual(testCase, cycling(array, 5, 5, true), [3, 2, 2, 1, 1]);
    
    end
    
    function divide_searching_set_test(testCase)
    %DIVIDE_SEARCHING_SET_TEST tests the file private/divide_searching_set.m
    m = 11;
    nb = 3;
    INDEX_SEARCHING_SET = cell(1,nb);
    INDEX_SEARCHING_SET{1} = [1, 2, 3, 4, 3, 2, 1];
    INDEX_SEARCHING_SET{2} = [5, 6, 7, 8, 7, 6, 5];
    INDEX_SEARCHING_SET{3} = [9, 10, 11, 10, 9];
    cycling = 5;
    polling = "opportunistic";
    
    verifyEqual(testCase, divide_searching_set(m, nb, cycling, polling), INDEX_SEARCHING_SET)
    
    end
    
    function get_default_constant_test(testCase)
    %GET_DEFAULT_CONSTANT_TEST tests the file private/get_default_constant.m
    constant_name = "expand";
    constant_value = 2;
    verifyEqual(testCase, get_default_constant(constant_name), constant_value)
    
    constant_name = "shrink";
    constant_value = 0.5;
    verifyEqual(testCase, get_default_constant(constant_name), constant_value)
    
    constant_name = "maxfun";
    constant_value = 1e4;
    verifyEqual(testCase, get_default_constant(constant_name), constant_value)
    
    constant_name = "sufficient_decrease_factor";
    constant_value = 1e-3;
    verifyEqual(testCase, get_default_constant(constant_name), constant_value)
    
    constant_name = "StepTolerance";
    constant_value = 1e-12;
    verifyEqual(testCase, get_default_constant(constant_name), constant_value)
    
    constant_name = "ftarget";
    constant_value = -inf;
    verifyEqual(testCase, get_default_constant(constant_name), constant_value)
    
    constant_name = "polling";
    constant_value = "opportunistic";
    verifyEqual(testCase, get_default_constant(constant_name), constant_value)
    
    constant_name = "polling_inner";
    constant_value = "opportunistic";
    verifyEqual(testCase, get_default_constant(constant_name), constant_value)
    
    constant_name = "polling_outer";
    constant_value = "opportunistic";
    verifyEqual(testCase, get_default_constant(constant_name), constant_value)
    
    constant_name = "cycling";
    constant_value = 3;
    verifyEqual(testCase, get_default_constant(constant_name), constant_value)
    
    constant_name = "memory";
    constant_value = true;
    verifyEqual(testCase, get_default_constant(constant_name), constant_value)
    
    end
    
    function get_exitflag_test(testCase)
    %GET_EXITFLAG_TEST tests the file private/get_exitflag.m
    information = "MAXIT_REACHED";
    EXITFLAG = 3;
    
    verifyEqual(testCase, get_exitflag(information), EXITFLAG)
    
    end
    
    % function inner_direct_search_test(testCase)
    % %INNER_DIRECT_SEARCH_TEST tests the file private/inner_direct_search.m
    % 
    % % TODO: Create the tests.
    % 
    % end
    
    function searching_set_test(testCase)
    %SEARCHING_SET_TEST tests the file private/searching_set.m
    n = 5;
    D = [zeros(n) zeros(n)];
    for i = 1:n
        D(i, 2*i-1) = 1;
        D(i, 2*i) = -1;
    end
    %D = [eye(n) -eye(n)];
    options.direction = "canonical";
    
    verifyEqual(testCase, searching_set(n,options), D)
    
    end
    
    % function blockwise_direct_search_test(testCase)
    % %BLOCKWISE_DIRECT_SEARCH_TEST tests the file private/blockwise_direct_search.m
    % 
    % end