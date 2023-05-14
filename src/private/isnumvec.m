function is_numvec = isnumvec(x)
% Verify whether x is a number vector.
is_numvec = isnumeric(x) && isvector(x);
end