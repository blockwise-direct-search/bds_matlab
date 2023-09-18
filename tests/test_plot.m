% 原始变量名
originalVariable = figure;

% 指定的新变量名
newVariableName = 'myNewVariable';

% 使用 eval 函数修改变量名
eval([newVariableName ' = originalVariable;']);

% 测试新变量名
disp(myNewVariable);