function verify_postconditions(fun, xopt, fopt, exitflag, output)
%VERIFY_POSTCONDITIONS verifies whether output is in right form.
%

% Verify whether xopt is a real vector.
if ~isrealvector(xopt)
    error("xopt is not a real vector.");
end

% Verify whether fopt is a real number.
if ~(isrealscalar(fopt))
    error("fopt is not a real number.");
end

% Verify whether exitflag is an integer.
if ~(isintegerscalar(exitflag))
    error("exitflag is not an integer.");
end

% Verify whether nf is a positive integer.
if ~isfield(output, "funcCount")
    error("output.funcCount does not exist.");
end
nf = output.funcCount;
if ~(isintegerscalar(nf) && nf > 0)
    error("output.funcCount is not a positive integer.");
end

% Verify whether output is a structure.
if ~(isstruct(output))
    error("output is not a structure.");
end

% Verify whether output.fhist exists.
if ~isfield(output, "fhist")
    error("output.fhist does not exist.");
end
fhist = output.fhist;
nhist = length(fhist);

if isfield(output, "xhist")
    xhist = output.xhist;
    % Verify whether xhist is a real matrix of size.
    if ~(isrealmatrix(xhist) && any(size(xhist) == [length(xopt), nhist]))
        error("output.xhist is not a real matrix.");
    end
    
    % Check whether length(fhist) is equal to length(xhist) and nf respectively.
    if ~(length(fhist) == size(xhist, 2) && size(xhist, 2) == nf)
        error("length of fhist is not equal to length of xhist or nf.");
    end
    
    % Check whether fhist == fun(xhist).
    fhist_eval = NaN(1, length(fhist));
    for i = 1:length(fhist)
        fhist_eval(i) = eval_fun(fun, xhist(:, i));
    end
    
    % In case of fhist_eval(i) = NaN or fhist(i) = NaN.
    % assert(all( (isnan(A) & isnan(B)) | A==B ));
    assert(all( (isnan(fhist) & isnan(fhist_eval)) | fhist==fhist_eval ));
    
    % Examine whether fopt is the minimum of fhist and xopt is the
    % corresponding point (when receiving simple decrease)
    % assert(fun(xopt) == fopt && min(fhist) == fopt);
end

end
