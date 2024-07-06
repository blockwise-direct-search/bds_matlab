% Defining the range of independent variables.
x1 = -10:0.1:10;
x2 = -10:0.1:10;
[X1, X2] = meshgrid(x1, x2);

% Calculate function value.
Y = quadratic([X1(:), X2(:)]);
Y = reshape(Y, size(X1));

% Draw a contour map.
figure;
contourf(X1, X2, Y, 20); % Draw 20 contour lines.
colorbar;
hold on

options.output_xhist = true;
options.Algorithm = "cbds";
[xopt, fopt, exitflag, output] = bds(@quadratic, [5+randn(1), 5+randn(1)], options);

% Plot the trajectory of the optimization process.
plot(output.xhist(1,:), output.xhist(2,:), '-o', 'LineWidth', 2);

% Add the sequential label to each point.
% for i = 1:size(output.xhist, 2)
for i = 1:50
    text(output.xhist(1,i), output.xhist(2,i), num2str(i), 'VerticalAlignment','bottom', 'HorizontalAlignment','right');
end

xlabel('X');
ylabel('Y');
title('Points Trajectory');


function y = quadratic(x)
    A = [sqrt(2)/2, sqrt(2)/2; -sqrt(2)/2, sqrt(2)/2];
    y = sum(x .* (x * A'), 2);
end