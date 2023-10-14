function [action, wrong_input] = parse_input(argin)
%This file is cited from https://github.com/libprima/prima/blob/main/matlab/setup_tools/parse_input.m,
%which is written by Zaikun Zhang.
%PARSE_INPUT parses the input to the setup script.

% Compilation options.
action_list = {'all', 'bds', 'uninstall'};
action = 'compile';
wrong_input = false;

% Start the parsing to set `input_string` and `options`.
input_string = 'ALL';  % Default value for `input_string`.
if length(argin) > 1
    fprintf('\nSetup accepts at most one inputs.\n\n');
    wrong_input = true;
elseif length(argin) == 1
    if ischarstr(argin{1})
        input_string = argin{1};
    else
        fprintf('\nThe input to setup should be a string and/or a structure.\n\n');
        wrong_input = true;
    end
end

% Cast input_string to a character array in case it is a MATLAB string.
input_string = lower(char(input_string));

% Parse `input_string` to set `action`.
if ismember(input_string, action_list)
    action = input_string;
else
    wrong_input = true;
end

return
