% cd '/home/lhtian97/Documents/bds/tests/private';
% locate_matcutest
% cd ..
% p = macup('akiva');
% obj = ScalarFunction(p);
% test_options.is_noisy = false;
% test_options.noise_type = "gaussian";
% test_options.is_abs_noise = false;
% test_options.noise_level = 1e-3;
% r = 2;
% options.maxfun = 1e3;
% options.Algorithm = "newuoa";
% nlopt(@(x)obj.fun(x,test_options.is_noisy,r,test_options), p.x0, options)
% obj
% 
% 要创建的数组大小（示例）
arraySize = [1000, 100];  % 行数 x 列数

% 获取当前 Java 虚拟机的运行时实例
rt = java.lang.Runtime.getRuntime;

% 获取最大可用内存大小
maxMemory = rt.maxMemory;

% 计算要创建的数组的总字节数
arrayBytes = prod(arraySize) * 8;  % 假设每个元素为 8 字节

% 检查数组大小是否超过最大数组大小限制
if arrayBytes <= maxMemory
    disp('数组大小在最大数组大小限制范围内。');
else
    disp('数组大小超过最大数组大小限制。');
end