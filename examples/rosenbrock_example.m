function rosenbrock_example()
%ROSENBROCK_EXAMPLE illustrates how to use prima.
%
%   ***********************************************************************
%   Authors:    LI Haitian (haitian-li@connect.polyu.hk) 
%               and Zaikun ZHANG (zaikun.zhang@polyu.edu.hk)
%               Department of Applied Mathematics,
%               The Hong Kong Polytechnic University
%
%   ***********************************************************************

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Attribute: public (can be called directly by users)
% 
% TODO: None
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\nMinimize the chained Rosenbrock function without constraints:\n');
x0 = [0; 0; 0];  % starting point

options.maxfun = 1e4;
options.StepTolerance = eps;
options.Algorithm = "cbds";

fullpath = mfilename('fullpath');
[path_examples,~] = fileparts(fullpath);
[path_bds, ~, ~] = fileparts(path_examples);
path_src = fullfile(path_bds, 'src');
addpath(path_src)

% The following syntax is identical to fmincon:
[xval, fval, exitflag, output] = bds(@rosenb, [0; 0; 0], options)
% Alternatively, the problem can be passed to bds as a structure:
%p.objective = @chrosen; p.x0 = x0;
%[x, fx, exitflag, output] = bds(p.objective, p.x0);

fprintf('\n1. No constraints:\n');
[xval, fval, exitflag, output] = bds(@chrosen, x0);

rmpath(path_src)

return


function f = chrosen(x)  % the subroutine defining the objective function
f = sum((x(1:end-1)-1).^2 + 4*(x(2:end)-x(1:end-1).^2).^2);
return
