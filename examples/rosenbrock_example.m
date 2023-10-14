function rosenbrock_example()
%ROSENBROCK_EXAMPLE illustrates how to use bds.
%
%   ***********************************************************************
%   Authors:    Haitian Li (haitian-li@connect.polyu.hk) 
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
path_examples = fileparts(fullpath);
path_bds = fileparts(path_examples);
path_src = fullfile(path_bds, 'src');
addpath(path_src)

fprintf('\n1. No constraints:\n');
% The following syntax is identical to fmincon:
[xval, fval, exitflag, output] = bds(@chrosen, x0, options);

rmpath(path_src)

return


function f = chrosen(x)  % the subroutine defining the objective function
f = sum((x(1:end-1)-1).^2 + 4*(x(2:end)-x(1:end-1).^2).^2);
return
