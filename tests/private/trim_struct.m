function parameters_saved = trim_struct(parameters_saved)
% Trim the string using "_" to replace " ", "-" and ":".

% 获取结构体的字段名
fields = fieldnames(parameters_saved);

% 将字段和对应的值写入文件
for i = 1:numel(fields)
    field = fields{i};
    value = parameters_saved.(field);
    value_length = length(value);
    if isnumvec(value)
        value = num2str(value);
        separator = ", ";
        if value_length ~= 1
            value = strjoin(strsplit(value), separator);
        end
        parameters_saved.(field) = value;
    elseif islogical(value) 
        if value
            parameters_saved.(field) = "true";
        else
            parameters_saved.(field) = "false";
        end
    elseif ischarstr(value)
        value = strjoin(value, ", ");
        parameters_saved.(field) = value;
    end
end

end
