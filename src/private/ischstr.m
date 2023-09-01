function is_chstr = ischstr(x)
%ISCHARSTR Check whether x is a 'char' or 'string'.

is_chstr = isa(x, "char") || isa(x, "string");