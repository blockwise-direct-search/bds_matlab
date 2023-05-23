function pmpaths = locate_prima(directory)
%This function finds where prima (https://github.com/libprima/prima.git) is installed, adds the
% paths needed for using prima, and returns these paths in a cell array.
% We search at most 3 levels below the given directory, whose default value is the home directory.

% We use the following path as the signature to identify prima.
signature_path = fullfile('prima', 'matlab', 'interfaces');

% cmtools is the path to the directory containing the signature path.
pmtools = '';

if nargin < 1
    path_strs = strsplit(path(), pathsep);
    ind = find(endsWith(path_strs, signature_path), 1, 'first');
    if ~isempty(ind)
        pmtools = path_strs{ind};
    else
        directory = getenv('HOME');
    end
end

if isempty(pmtools)
    % In the following line, the "*/" before signature_path cannot be removed.
    name_str = ['"*/', signature_path, '"'];
    [~, pmtools] = system(['find ', directory, ' -maxdepth 9 -wholename ', name_str, ' -type d -print -quit']);

    if isempty(pmtools)
        error('locate_prima:primaNotFound', 'Prima is not found under %s.', directory);
    end
end

pmtools = strtrim(pmtools);  % Remove the leading and trailing white-space characters, including '\n'.

pmpaths = {pmtools};  % There may be other paths to include in the future.

for ip = 1 : length(pmpaths)
    addpath(pmpaths{ip});
end
