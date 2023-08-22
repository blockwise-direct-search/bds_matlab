function iscs = ischarstr(x)
%ISCHARSTR Check whether an input is a 'char' or 'string'.

iscs = (isa(x, "char") || isa(x, "string"));