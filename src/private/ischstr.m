function is_chstr = ischstr(x)
% This file is cited from https://github.com/libprima/prima/blob/main/matlab/tests/private/ischstr.m, which is
% written by Zaikun ZHANG.
%ISCHARSTR checks whether x is a "char" or "string".

is_chstr = isa(x, "char") || isa(x, "string");