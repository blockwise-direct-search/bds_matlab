cd '/home/lhtian97/Documents/bds/tests/private';
locate_matcutest
cd ..
p = macup('ALLINITU');
obj = ScalarFunction(p);
test_options.is_noisy = false;
test_options.noise_type = "gaussian";
test_options.is_abs_noise = false;
test_options.noise_level = 1e-3;
r = 2;
options.maxfun = 1e5;
fminsearch(@(x)obj.fun(x,test_options.is_noisy,r,test_options), p.x0, options)