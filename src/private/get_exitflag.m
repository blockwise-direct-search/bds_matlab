function [exitflag] = get_exitflag(information)
%GET_EXITFLAG gets the EXITFLAG of BDS.
%   SMALL_ALPHA     Step size is below StepTolerance. In the case of variable step sizes, 
%                   SMALL_ALPHA indicates the largest component of step sizes is below StepTolerance.
%   MAXFUN_REACHED  The number of function evaluations reaches MAXFUN.
%   FTARGET_REACHED Function value is smaller than or equal to FTARGET.
%   MAXIT_REACHED   The number of iterations reaches MAXIT.  
%

% Check whether INFORMATION is a string or not.
if is_debugging
    if ~isstring(information)
        error("Information is not a string.");
    end
end

switch information
    case "SMALL_ALPHA"
        exitflag = 0;
    case "FTARGET_REACHED"
        exitflag = 1;
    case "MAXFUN_REACHED"
        exitflag = 2;
    case "MAXIT_REACHED"
        exitflag = 3;
    otherwise
        exitflag = -1;
end

%exitflag = find(break_conditions == information) - 1;
if isempty(exitflag)
    exitflag = -1;
    disp("New break condition happens."); 
end

% Check whether EXITFLAG is an integer or not.
if is_debugging
    if ~isintegerscalar(exitflag)
        error("Exitflag is not an integer.");
    end
end 

end
