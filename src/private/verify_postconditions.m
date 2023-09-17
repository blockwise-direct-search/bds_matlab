function verify_postconditions(fun, xval, fval, exitflag, output)
%VERIFY_POSTCONDITIONS verifies whether output is in right form.
%

% Verify whether xval is a real column or scalar.
if ~(isnumeric(xval) && isreal(xval) && iscolumn(xval))
    error("xval is not a real column or scalar.");
end

% Verify whether fval is a real number.
if ~(isrealscalar(fval))
    error("fval is not a real number.");
end

% Verify whether exitflag is an integer.
if ~(isintegerscalar(exitflag))
    error("exitflag is not an integer.");
end

% Verify whether nf is a postive integer.
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
    if ~(isrealmatrix(xhist) && any(size(xhist) == [length(xval), nhist]))
        error("output.xhist is not a real matrix.");
    end
    
    % Check whether length(fhist) is equal to length(xhist) and nf respectively.
    if ~(length(fhist) == size(xhist, 2) && size(xhist, 2) == nf)
        error("length of fhist is not equal to length of xhist or nf.");
    end
    
    % Check whether fhist == fun(xhist).
    % TODO: there is a way to avoid using loop, try to find it.
    fhist_eval = NaN(1, length(fhist));
    for i = 1:length(fhist)
        fhist_eval(i) = eval_fun(fun, xhist(:, i));
    end
    
    % In case of fhist_eval(i) = NaN or fhist(i) = NaN.
    % assert(all( (isnan(A) & isnan(B)) | A==B ));
    assert(all( (isnan(fhist) & isnan(fhist_eval)) | fhist==fhist_eval ));
    
    % Examine whether fval is the minimum of fhist and xval is the
    % corresponding point (when receiving simple decrease)
    % assert(fun(xval) == fval && min(fhist) == fval);
end

end
