restoredefaultpath;
% The code of the following lines is for using matcutest first time.
%addpath('/home/htl/local/matcutest/mtools/src');
addpath('/home/lhtian97/local/matcutest/mtools/src');
addpath('/home/lhtian97/bds_new_framework/tests/competitors/prima/matlab/interfaces/');
fullpath = mfilename('fullpath');
[path_tests,~] = fileparts(fullpath);
% The code of the following three lines is for running prima first time.
% cd(path_tests)
% cd ./competitors/prima
% setup
addpath(path_tests);
cd(path_tests)
cd ..
path_bds = pwd;
addpath(path_bds);
path_src = strcat(path_bds, "/src");
addpath(path_src);
path_competitors = strcat(path_tests, "/competitors");
addpath(path_competitors);
parameters.solvers_invoke = ["blockwise_direct_search",  ...
  "newuoa"];
parameters.solvers_label = [ ...
 "GS", "newuoa"];
parameters.memory = [true, true];
parameters.polling_outer = ["opportunistic", "opportunistic",...
    ];
parameters.polling_blocks = ["Gauss-Seidel", "Randomized_no_repetition"];
parameters.polling_inner = ["opportunistic", "opportunistic",...
    ];
parameters.cycling_inner = [1, 1];
parameters.solvers_tag = ["GS", "newuoa"];
parameters.nb_generator = [0.5, 0.5];
nb_tag = ["n", "n"];
parameters.maxfun_dim = 1000;
parameters.maxfun = 60000;
parameters.problems_type = 'u';
parameters.problems_mindim = 6;
parameters.problems_maxdim = 60;
parameters.tau = 10.^(-1:-1:-10);
% canonical; random
parameters.direction = ["canonical", "canonical"];
parameters.parallel = false;
num_solvers = length(parameters.cycling_inner);
pdfname = "";
% Name pdf automatically (not manually).
for i = 1:num_solvers
    if i > 1
       pdfname = strcat(pdfname, "_"); 
    end
    pdfname = strcat(pdfname, parameters.solvers_tag(i), "_", nb_tag(i),...
         "_",  num2str(parameters.cycling_inner(i)));
end
pdfname = strcat(pdfname, "_", num2str(parameters.problems_mindim), "_", num2str(parameters.problems_maxdim));
parameters.pdfname = pdfname;
testbds(parameters);
% Delete the path to recover
rmpath(path_tests);
rmpath(path_bds);
rmpath(path_src);
rmpath(path_competitors);

