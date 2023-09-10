function wm_newuoa(fun,x0,options)
% This is a wrapper for NEWUOA algorithm which was implemented by Stefan
% Wild and Jorge More.
%

currentFile = mfilename('fullpath');
currentPath = fileparts(currentFile);
path_mnewuoa = fullfile(currentPath, 'mnewuoa');
addpath(path_mnewuoa);

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

mnewuoa(fun,n,options.npt,x0',...
    options.rhobeg,options.rhoend,options.iprint,options.maxfun);

end

