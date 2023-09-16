function tests = unit_test
    %   UNIT_TEST runs all the test functions in this file.
    %   To run these tests, simply type "Run Tests" in the command window. To
    %   create a new test function in this file with a name that starts or
    %   finishes with "test" (case insensitive). For more info, see
    %
    %   https://www.mathworks.com/help/matlab/matlab_prog/write-simple-test-case-with-functions.html
    
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
    
    % The following are the ones commented on cycling.m.
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
    
    end
    
    function divide_searching_set_test(testCase)
    %DIVIDE_SEARCHING_SET_TEST tests the file private/divide_searching_set.m
    m = 11;
    nb = 3;
    INDEX_SEARCHING_SET = cell(1,nb);
    INDEX_SEARCHING_SET{1} = [1, 2, 3, 4];
    INDEX_SEARCHING_SET{2} = [5, 6, 7, 8];
    INDEX_SEARCHING_SET{3} = [9, 10, 11];
    
    verifyEqual(testCase, divide_searching_set(m, nb), INDEX_SEARCHING_SET)
    
    end
    
    function get_default_constant_test(testCase)
    %GET_DEFAULT_CONSTANT_TEST tests the file private/get_default_constant.m
    constant_name = "maxfun";
    constant_value = 1e5;
    verifyEqual(testCase, get_default_constant(constant_name), constant_value)

    constant_name = "maxfun_factor";
    constant_value = 1e3;
    verifyEqual(testCase, get_default_constant(constant_name), constant_value)

    constant_name = "expand";
    constant_value = 2;
    verifyEqual(testCase, get_default_constant(constant_name), constant_value)
    
    constant_name = "shrink";
    constant_value = 0.5;
    verifyEqual(testCase, get_default_constant(constant_name), constant_value)
    
    constant_name = "sufficient_decrease_factor";
    constant_value = 1e-3;
    verifyEqual(testCase, get_default_constant(constant_name), constant_value)

    constant_name = "accept_simple_decrease";
    constant_value = false;
    verifyEqual(testCase, get_default_constant(constant_name), constant_value)    
     
    constant_name = "StepTolerance";
    constant_value = eps;
    verifyEqual(testCase, get_default_constant(constant_name), constant_value)

    constant_name = "shuffle_period";
    constant_value = 1;
    verifyEqual(testCase, get_default_constant(constant_name), constant_value)

    constant_name = "replacement_delay";
    constant_value = 0;
    verifyEqual(testCase, get_default_constant(constant_name), constant_value)
     
    constant_name = "ftarget";
    constant_value = -inf;
    verifyEqual(testCase, get_default_constant(constant_name), constant_value)
     
    constant_name = "polling_inner";
    constant_value = "opportunistic";
    verifyEqual(testCase, get_default_constant(constant_name), constant_value)
    
    constant_name = "with_cycling_memory";
    constant_value = true;
    verifyEqual(testCase, get_default_constant(constant_name), constant_value)
    
    end
    
    function get_exitflag_test(testCase)
    %GET_EXITFLAG_TEST tests the file private/get_exitflag.m.

    information = "SMALL_ALPHA";
    EXITFLAG = 0;    
    verifyEqual(testCase, get_exitflag(information), EXITFLAG)

    information = "MAXFUN_REACHED";
    EXITFLAG = 1;    
    verifyEqual(testCase, get_exitflag(information), EXITFLAG)

    information = "FTARGET_REACHED";
    EXITFLAG = 2;    
    verifyEqual(testCase, get_exitflag(information), EXITFLAG)

    information = "MAXIT_REACHED";
    EXITFLAG = 3;    
    verifyEqual(testCase, get_exitflag(information), EXITFLAG)
    
    end
    
    function inner_direct_search_test(testCase)
    %INNER_DIRECT_SEARCH_TEST tests the file private/inner_direct_search.m.
    
    fun = @(x) x(1)^2 + x(2)^2;
    xval = [1; 1];
    fval = 2;
    D = [1 -1; 0 0];
    direction_indices = [1 2];
    alpha = 1;
    options.maxfun = 2;

    xval_result = [0; 1];
    fval_result = 1;
    exitflag_result = NaN;
    direction_indices_result = [2 1];
    fhist_result = [5 1];
    xhist_result = [2 0;1 1];
    nf_result = 2;
    success_result = true;
    terminate_result = false;

    [xval_update, fval_update, exitflag, output] = inner_direct_search(fun, ...
    xval, fval, D, direction_indices, alpha, options);
    
    verifyEqual(testCase, xval_update, xval_result);
    verifyEqual(testCase, fval_update, fval_result);
    verifyEqual(testCase, exitflag, exitflag_result);

    verifyEqual(testCase, output.fhist, fhist_result);
    verifyEqual(testCase, output.xhist, xhist_result);
    verifyEqual(testCase, output.nf, nf_result);
    verifyEqual(testCase, output.success, success_result);
    verifyEqual(testCase, output.direction_indices, direction_indices_result);
    verifyEqual(testCase, output.terminate, terminate_result);
    

    end
    
    function searching_set_test(testCase)
    %SEARCHING_SET_TEST tests the file private/searching_set.m.
    n = 5;
    D = [zeros(n) zeros(n)];
    for i = 1:n
        D(i, 2*i-1) = 1;
        D(i, 2*i) = -1;
    end
    options.direction = "canonical";    
    verifyEqual(testCase, get_searching_set(n,options), D)

    options = struct();
    verifyEqual(testCase, get_searching_set(n, options), D)
    
    options.direction = "identity";
    D = [eye(n) -eye(n)];
    verifyEqual(testCase, get_searching_set(n, options), D)



    end
    
    function blockwise_direct_search_test(testCase)
    %BLOCKWISE_DIRECT_SEARCH_TEST tests the file bds.m.

    fun = @(x) x(1)^2 + x(2)^2;
    x0 = [100; 100];
    options = struct();

    xval_result = [0; 0];
    fval_result = 0;
    exitflag_result = 0;

    [xval, fval, exitflag] = bds(fun, x0, options);

    verifyEqual(testCase, xval, xval_result);
    verifyEqual(testCase, fval, fval_result);
    verifyEqual(testCase, exitflag, exitflag_result);
    
    end