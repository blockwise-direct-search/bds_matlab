function [outputArg1,outputArg2] = prof(inputArg1,inputArg2)
% 生成一些类型为figure的变量
figure1 = figure;
plot([1, 2, 3], [4, 5, 6]);
title('Figure 1');

figure2 = figure;
bar([10, 20, 30]);
title('Figure 2');

figure3 = figure;
plot([1, 4, 9], [2, 5, 10]);
title('Figure 3');

% 创建一个新的Figure用于合并
mergedFigure = figure;

% 定义要合并的figure变量
figures = [figure1, figure2, figure3];

% 获取要合并的figure的数量
numFigures = numel(figures);

% 设置子图的行数和列数
numRows = 1;
numCols = numFigures;

% 循环处理每个figure
for i = 1:numFigures
    % 将当前figure放置在对应的子图中
    subplot(numRows, numCols, i);
    copyobj(allchild(get(figures(i), 'CurrentAxes')), gca);
end

% 显示合并后的Figure
title('Merged Figure');

% 将合并后的figure保存为EPS格式
saveas(mergedFigure, 'merged_figure.eps', 'epsc');
end

