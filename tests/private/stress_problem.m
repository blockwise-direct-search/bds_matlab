function problem = stress_problem(n, problem_type, random_seed)
    %This function generates a test problem for the stress test.
    % n: the dimension of the problem
    % problem_type:
    %   'u' for unconstrained,
    % random_seed: the random seed for the problem generation
    
    % Set the random seed
    orig_rng_state = rng();
    rng(random_seed);
    
    % Set the starting point
    problem.x0 = randn(n, 1);
    problem.objective = @chrosen;
    
    % Restore the random seed
    rng(orig_rng_state);
    
    % `stress_problem` ends here
    return   
    
    function f = chrosen(x)  % the subroutine defining the objective function
    alpha = 4;
    f = sum((x(1:end-1)-1).^2 + alpha*(x(2:end)-x(1:end-1).^2).^2);
    return