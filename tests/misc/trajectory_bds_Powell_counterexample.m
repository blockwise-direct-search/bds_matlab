% value = Powell_counterexample_1([-1, 1, -1])
% keyboard
options.MaxFunctionEvaluations = 1e6;
options.output_xhist = true;
options.Algorithm = "cbds";
[xopt, fopt, exitflag, output] = bds(@Powell_counterexample_1, [-1-1e-6, 1+0.5*1e-6, -1-0.25*1e-6], options);
% fopt == Powell_counterexample_1(xopt)


% 创建示例矩阵
points = output.xhist;

% 提取点的坐标
x = points(1, :);
y = points(2, :);
z = points(3, :);

% 绘制散点图
scatter3(points(1, :), points(2, :), points(3, :), 'filled');

% 设置图形属性
xlabel('X');
ylabel('Y');
zlabel('Z');
title('3D Scatter Plot');
grid on;

function y = Powell_counterexample_1(x)

y = -x(1)*x(2) - x(2)*x(3) - x(3)*x(1) + Powell_counterexample_basis(x(1), 1) ...
    + Powell_counterexample_basis(-x(1), 1) + Powell_counterexample_basis(x(2), 1) ...
    + Powell_counterexample_basis(-x(2), 1) + ...
    + Powell_counterexample_basis(x(3), 1) + Powell_counterexample_basis(-x(3), 1);

end

function y = Powell_counterexample_basis(x, c)

if x <= c
    y = 0;
else
    y = (x - c)^2;
end

end