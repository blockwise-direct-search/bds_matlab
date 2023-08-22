function [exitflag] = get_exitflag(information)
%To get exitflag of variable situations.
%     SMALL_ALPHA     Step size is below StepTolerance. For the case of variable
%                     step sizes, it indicates the maximum of step sizes is below
%                     StepTolerance.
%     MAXFUN_REACHED  The number of function evaluations equal to the maxfun.
%     FTARGET_REACHED Function value is less or equal to ftarget.
%     MAXIT_REACHED   The number of iterations equal to maxit.  

% Preconditions: information should be a string.
if is_debugging
    assert(isstring(information));
end

break_conditions = ["SMALL_ALPHA";"MAXFUN_REACHED";"FTARGET_REACHED";"MAXIT_REACHED"];

exitflag = find(break_conditions == information) - 1;
if isempty(exitflag)
    exitflag = -1;
    disp("New break condition happens."); 
end

% Postcondtions: exitflag should be an integer.
if is_debugging
    assert(ceil(exitflag) == exitflag);
end 

end
