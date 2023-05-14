function [result] = trim_time(str)
% Trim the format of string use '_'
         result = replace(str,' ','_');
         result = replace(result,'-','_');
         result = replace(result,':','_');
end
