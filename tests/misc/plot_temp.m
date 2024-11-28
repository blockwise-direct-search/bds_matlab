function plot_temp(x, y, z)
    % plotBlueGradientSurface 绘制蓝色渐变插值表面图
    % 输入:
    %   x - x 坐标的向量，范围应为 [1, 10]
    %   y - y 坐标的向量，范围应为 [0, 1]
    %   z - z 坐标的三维数组，大小为 [length(x), length(y), n]
    
    % 创建均匀网格
    [xq, yq] = meshgrid(linspace(1, 10, 100), linspace(0, 1, 100)); % 生成均匀网格

    % 将 x 和 y 的网格点展开为向量
    [X, Y] = meshgrid(x, y); % 创建 x 和 y 的网格
    X = X(:); % 将 x 网格转换为列向量
    Y = Y(:); % 将 y 网格转换为列向量

    % 将 z 数据重塑为一维
    Z = reshape(z, numel(X), []); % 将 z 变为每组数据的形式

    % 使用插值计算对应的 z 值
    zq = nan(size(xq, 1), size(xq, 2), size(Z, 2)); % 初始化插值结果
    
    for i = 1:size(Z, 2)
        zq(:,:,i) = griddata(X, Y, Z(:, i), xq, yq, 'cubic'); % 对每组 z 进行插值
    end
    keyboard
    % 绘制插值后的图形
    figure; % 创建新图形窗口
    hold on; % 保持当前图形

    % 绘制每组插值后的表面图
    for i = 1:size(zq, 3)
        surf(xq, yq, zq(:,:,i), 'EdgeColor', 'none', 'FaceColor', 'interp'); % 绘制插值后的表面图
    end

    % 自定义颜色映射为蓝色渐变
    cmap = [linspace(0, 0, 256)', linspace(0, 0, 256)', linspace(0, 1, 256)']; % 生成蓝色渐变
    colormap(cmap); % 应用自定义颜色映射

    % 设置颜色条
    colorbar; % 显示颜色条
    xlabel('X-axis'); % x 轴标签
    ylabel('Y-axis'); % y 轴标签
    zlabel('Z-axis'); % z 轴标签
    title('Interpolated Surfaces with Blue Color Gradient'); % 图形标题
    view(2); % 从上方查看
    grid on; % 显示网格
end
