function tests = unit_test
    %   UNIT_TEST runs all the test functions in this file.
    %   To run these tests, simply type "Run Tests" in the command window. 
    %   To create a new test function in this file with a name that starts or
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
    
    function divide_direction_set_test(testCase)
    %DIVIDE_direction_set_TEST tests the file private/divide_direction_set.m
    n = 11;
    nb = 3;
    INDEX_direction_set = cell(1,nb);
    INDEX_direction_set{1} = [1 2 3 4 5 6 7 8];
    INDEX_direction_set{2} = [9 10 11 12 13 14 15 16];
    INDEX_direction_set{3} = [17 18 19 20 21 22];

    verifyEqual(testCase, divide_direction_set(n, nb), INDEX_direction_set)

    n = 10;
    nb = 3;
    INDEX_direction_set = cell(1,nb);
    INDEX_direction_set{1} = [1 2 3 4 5 6 7 8];
    INDEX_direction_set{2} = [9 10 11 12 13 14];
    INDEX_direction_set{3} = [15 16 17 18 19 20];
    
    verifyEqual(testCase, divide_direction_set(n, nb), INDEX_direction_set)

    n = 15;
    nb = 3;
    INDEX_direction_set = cell(1,nb);
    INDEX_direction_set{1} = [1 2 3 4 5 6 7 8 9 10];
    INDEX_direction_set{2} = [11 12 13 14 15 16 17 18 19 20];
    INDEX_direction_set{3} = [21 22 23 24 25 26 27 28 29 30];
    
    verifyEqual(testCase, divide_direction_set(n, nb), INDEX_direction_set)

    n = 3;
    nb = 3;
    INDEX_direction_set = cell(1,nb);
    INDEX_direction_set{1} = [1 2];
    INDEX_direction_set{2} = [3 4];
    INDEX_direction_set{3} = [5 6];

    verifyEqual(testCase, divide_direction_set(n, nb), INDEX_direction_set)
    
    end
    
    function get_default_constant_test(testCase)
    %GET_DEFAULT_CONSTANT_TEST tests the file private/get_default_constant.m.
    constant_name = "MaxFunctionEvaluations_dim_factor";
    constant_value = 500;
    verifyEqual(testCase, get_default_constant(constant_name), constant_value)

    constant_name = "expand";
    constant_value = 2;
    verifyEqual(testCase, get_default_constant(constant_name), constant_value)
    
    constant_name = "shrink";
    constant_value = 0.5;
    verifyEqual(testCase, get_default_constant(constant_name), constant_value)

    assert(strcmp(func2str(get_default_constant("forcing_function")), func2str(@(alpha) alpha^2)));

    constant_name = "reduction_factor";
    constant_value = [0, eps, eps];
    verifyEqual(testCase, get_default_constant(constant_name), constant_value)  
     
    constant_name = "StepTolerance";
    constant_value = 1e-6;
    verifyEqual(testCase, get_default_constant(constant_name), constant_value)
    
    constant_name = "permuting_period";
    constant_value = 1;
    verifyEqual(testCase, get_default_constant(constant_name), constant_value)

    constant_name = "replacement_delay";
    constant_value = 1;
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

    information = "FTARGET_REACHED";
    EXITFLAG = 1;    
    verifyEqual(testCase, get_exitflag(information), EXITFLAG)

    information = "MAXFUN_REACHED";
    EXITFLAG = 2;    
    verifyEqual(testCase, get_exitflag(information), EXITFLAG)

    information = "MAXIT_REACHED";
    EXITFLAG = 3;    
    verifyEqual(testCase, get_exitflag(information), EXITFLAG)
    
    end
    
    function direction_set_test(testCase)
    %direction_set_TEST tests the file private/direction_set.m.
    n = 5;
    D = [zeros(n) zeros(n)];
    for i = 1:n
        D(i, 2*i-1) = 1;
        D(i, 2*i) = -1;
    end
    verifyEqual(testCase, get_direction_set(n), D)

    options = struct();
    verifyEqual(testCase, get_direction_set(n, options), D)

    n = 3;
    options = struct();
    A = randn(n);
    [Q, ~] = qr(A);
    options.direction_set = Q;
    D = get_direction_set(n, options);
    if D(:, 1:2:5) ~= -D(:, 2:2:6)
        error('D is not symmetric');
    end    
    if rank(D(:, 1:2:5)) ~= 3
        error('The basis of D is not a basis.');
    end

    n = 3;
    options = struct();
    A = randn(n);
    [Q, ~] = qr(A);
    options.direction_set = Q;
    D = get_direction_set(n, options);
    if D(:, 1:2:5) ~= -D(:, 2:2:6)
        error('D is not symmetric');
    end    
    if rank(D(:, 1:2:5)) ~= 3
        error('The basis of D is not a basis.');
    end

    n = 3;
    options = struct();
    A = zeros(3);
    detA = 0;
    while detA == 0
            A = randn(3);
            detA = det(A);
    end
    options.direction_set = A;
    D = get_direction_set(n, options);
    if D(:, 1:2:5) ~= -D(:, 2:2:6)
        error('D is not symmetric');
    end    
    if rank(D(:, 1:2:5)) ~= 3
        error('The basis of D is not a basis.');
    end
  
    end
    