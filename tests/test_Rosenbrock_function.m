function test_Rosenbrock_function()
%This file is cited from https://github.com/libprima/prima/blob/main/matlab/examples/rosenbrock_example.m, which is
%written by Zaikun Zhang.
%ROSENBROCK_EXAMPLE illustrates how to use bds.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Attribute: public (can be called directly by users)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%rng(267);
%x0 = [0; 0; 0];  % starting point
x0 = randn(20,1);  % starting point

%options.MaxFunctionEvaluations = 1e4;
%options.StepTolerance = eps;
options = struct();
options.StepTolerance = eps;

fullpath = mfilename("fullpath");
path_examples = fileparts(fullpath);
path_bds = fileparts(path_examples);
path_src = fullfile(path_bds, "src");
path_competitors = fullfile(path_bds, "tests", "competitors");
addpath(path_src)
addpath(path_competitors)

% The following syntax is identical to fmincon:
%[X,FVAL,EXITFLAG,OUTPUT] = fminsearch(@chrosen, x0)
%[X,FVAL,EXITFLAG,OUTPUT] = fminunc(@chrosen, x0)
[xopt, fopt, exitflag, output] = bds(@chrosen, x0)
bfo_wrapper(@chrosen, x0, options)
rmpath(path_src)
rmpath(path_competitors)
return

function f = chrosen(x)  % the subroutine defining the objective function
f = sum((x(1:end-1)-1).^2 + 4*(x(2:end)-x(1:end-1).^2).^2);
% f = f*(1+1e-8*randn(1));
return

