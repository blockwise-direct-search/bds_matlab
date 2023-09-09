function [xval, fval, exitflag, output] = wm_newuoa(fun,x0,options)
% This is a wrapper for NEWUOA algorithm which was implemented by Stefan
% Wild.
%

currentFile = mfilename('fullpath');
[currentPath, ~, ~] = fileparts(currentFile);
path_mnewuoa = fullfile(currentPath, 'mnewuoa');
addpath(path_mnewuoa);

global mnewuoa_fhist
global mnewuoa_nf

% Dimension
n = numel(x0);

% Set options to an empty structure if it is not supplied.
if nargin < 3
    options = struct();
end

if ~isfield(options, 'npt')
    options.npt = 2*n+1;
end

if ~isfield(options, 'rhobeg')
    options.rhobeg = 1;
end

if ~isfield(options, 'rhoend')
    options.rhoend = 1e-6;
end

if ~isfield(options, 'iprint')
    options.iprint = 0;
end

if ~isfield(options, 'maxfun')
    options.maxfun = 1000*n;
end

mnewuoa_nf = 0;
mnewuoa_fhist = NaN(1, options.maxfun);

xval = mnewuoa(@(x)objective(fun,x'),n,options.npt,x0',...
    options.rhobeg,options.rhoend,options.iprint,options.maxfun);
xval = xval';

fval = fun(xval);

exitflag = -1;

output.nf = mnewuoa_nf;
output.fhist = mnewuoa_fhist(1:mnewuoa_nf);

end

function [f] = objective(fun, x)

global mnewuoa_fhist
global mnewuoa_nf
f = fun(x);
mnewuoa_nf = mnewuoa_nf+1;
mnewuoa_fhist(mnewuoa_nf) = f;

end
