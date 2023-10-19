% 创建一个示例表格
data = {'John', 25, 180.5;
        'Alice', 30, 165.2;
        'Bob', 35, 175.0};
    
columnNames = {'Name', 'dimension', 'ratio'};
    
T = cell2table(data, 'VariableNames', columnNames);

% 将表格保存到文本文件
writetable(T, 'table_file.txt', 'Delimiter', '\t');