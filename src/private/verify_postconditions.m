function verify_postconditions(fun, xval, fval, exitflag, output)

% Verify whether xval is a real column or scalar
assert(isnumeric(xval) && isreal(xval) && iscolumn(xval));

% Verify whether fval is a real number
assert(isrealscalar(fval));

% Verify whether exitflag is an integer
% TODO: examine whether the value of exitflag is corresponding to the
% information
assert(isintegerscalar(exitflag));

% Verify whether nf is a postive integer
assert(isfield(output, "funcCount"));
nf = output.funcCount;
assert(isintegerscalar(nf) && nf > 0);

% Verify whether output is a structure
assert(isa(output, "struct"));

% Verify whether output.fhist exists
assert(isfield(output, "fhist"));
fhist = output.fhist;
nhist = length(fhist);

% Verify whether output.xhist exists
assert(isfield(output, "xhist"));
xhist = output.xhist;
% Verify whether xhist is a real matrix of size
assert(isrealmatrix(xhist) && any(size(xhist) == [length(xval), nhist]));

% Check whether length(fhist) is equal to length(xhist) and nf
assert(size(xhist, 2) == nf && length(fhist) == nf);

% Check whether fhist = fun(xhist)
% TODO: there is a way to avoid using loop, try to find it.
fhistx = NaN(1, length(fhist));
for i = 1:length(fhist)
    fhistx(i) = fun(xhist(:, i));
end

% In case of fhistx(i) = NaN or fhist(i) = NaN.
% assert(all( (isana(A) & isnan(B)) | A==B ));
assert(all( (isnan(fhist) & isnan(fhistx)) | fhist==fhistx ));

% Examine whether fval is the minimum of fhist and xval is the
% corresponding point (when receiving simple decrease)
% assert(fun(xval) == fval && min(fhist) == fval);



% for i = 1:length(fhist)
%     fhistx(i) = fun(xhist(:, i));
%     % In case of fhistx(i) = NaN or fhist(i) = NaN.
%     assert(fhistx(i) == fhist(i) || (isnan(fhistx(i)) && isnan(fhist(i))));
%     if isnan(fhistx(i))
%         assert(isnan(isnan(fhist(i))));
%     else
%        assert(fhistx(i) == fhist(i));
%     end
%
% end

end
