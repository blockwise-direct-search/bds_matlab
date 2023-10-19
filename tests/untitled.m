% 创建一个示例表格
data = {'MUONSINELS', 1, 4.4983e+05, 0.0044983};
columnNames = {'Name', 'dimension', 'ratio', 'gval'};
T = cell2table(data, 'VariableNames', columnNames);

% 指定保存路径和文件名
filePath = './file.txt';

% 打开文本文件进行写入
fileID = fopen(filePath, 'w');

if fileID == -1
    error('Cannot open the file for writing.');
end

% 写入列名到文件并左对齐
fprintf(fileID, '%-15s\t%-15s\t%-15s\t%-15s\n', columnNames{:});

% 写入表格数据到文件并左对齐
for row = 1:size(T, 1)
    fprintf(fileID, '%-15s\t%-15d\t%-15g\t%-15g\n', T.Name{row}, T.dimension(row), T.ratio(row), T.gval(row));
end

% 关闭文件
fclose(fileID);