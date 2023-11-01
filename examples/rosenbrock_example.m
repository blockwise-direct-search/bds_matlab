function rosenbrock_example()
%This file is cited from https://github.com/libprima/prima/blob/main/matlab/examples/rosenbrock_example.m, which is
%written by Zaikun Zhang. 
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x0 = [0; 0; 0];  % starting point

% options.maxfun = 1e4;
% options.StepTolerance = eps;
% options.Algorithm = "cbds";
options.Algorithm = "newuoa";

fullpath = mfilename("fullpath");
path_examples = fileparts(fullpath);
path_bds = fileparts(path_examples);
path_src = fullfile(path_bds, "src");
path_competitors = fullfile(path_bds, "tests", "competitors");
addpath(path_src)
addpath(path_competitors)

% The following syntax is identical to fmincon:
% [xval, fval, exitflag, output] = nlopt_wrapper(@chrosen, x0, options)
[xopt, fmin, retcode] = nlopt_wrapper(@chrosen, x0, options)

rmpath(path_src)
rmpath(path_competitors)

return


function f = chrosen(x)  % the subroutine defining the objective function
f = sum((x(1:end-1)-1).^2 + 4*(x(2:end)-x(1:end-1).^2).^2);
return
