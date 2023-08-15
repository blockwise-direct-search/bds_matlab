function [f] = objective(x)
%OBJECTIVE Summary of this function goes here
%   Detailed explanation goes here
f = x' * x;
f = tough(f, x, 0);
end

