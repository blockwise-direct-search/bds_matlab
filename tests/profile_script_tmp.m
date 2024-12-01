clear all

parameters.parallel = false;

parameters.cut_margin = false;

% parameters.solvers_name = ["cbds", "newuoa"];

% parameters.problem_dim = "big";

% parameters.test_type = "s2mpj";

% parameters.multi_feature = ["plain", "noise_1e-1_no_rotation", "noise_1e-2_no_rotation", ...
%         "noise_1e-3_no_rotation", "noise_1e-4_no_rotation", "noise_1e-5_no_rotation", ...
%         "noise_1e-6_no_rotation", "noise_1e-7_no_rotation", "noise_1e-8_no_rotation"];

% profile_multi_feature(parameters);

% parameters.multi_feature = ["rotation", "rotation_noisy_1e-1",...
%         "rotation_noisy_1e-2", "rotation_noisy_1e-3", "rotation_noisy_1e-4",...
%         "rotation_noisy_1e-5", "rotation_noisy_1e-6", "rotation_noisy_1e-7",...
%         "rotation_noisy_1e-8"];

% profile_multi_feature(parameters);

% parameters.test_type = "matcutest";

% parameters.multi_feature = ["plain", "noise_1e-1_no_rotation", "noise_1e-2_no_rotation", ...
%         "noise_1e-3_no_rotation", "noise_1e-4_no_rotation", "noise_1e-5_no_rotation", ...
%         "noise_1e-6_no_rotation", "noise_1e-7_no_rotation", "noise_1e-8_no_rotation"];

% profile_multi_feature(parameters);

% parameters.multi_feature = ["rotation", "rotation_noisy_1e-1",...
%         "rotation_noisy_1e-2", "rotation_noisy_1e-3", "rotation_noisy_1e-4",...
%         "rotation_noisy_1e-5", "rotation_noisy_1e-6", "rotation_noisy_1e-7",...
%         "rotation_noisy_1e-8"];

% profile_multi_feature(parameters);


parameters.solvers_name = ["cbds", "nomad"];

parameters.problem_mindim = 2;
parameters.problem_maxdim = 5;
parameters.solvers_options{1}.expand = 1.25;
parameters.solvers_options{1}.shrink = 0.85;
% parameters.solvers_options{2}.expand = 2;
% parameters.solvers_options{2}.shrink = 0.5;

parameters.num_random = 1;

parameters.feature = "plain";

parameters.test_type = "s2mpj";

plot_profile(parameters);

% parameters.feature = "rotation_noisy_1e-3";

% plot_profile(parameters);

% parameters.feature = "noise_1e-3_no_rotation";

% plot_profile(parameters);

% parameters.multi_feature = ["randomx0_1e-3", "randomx0_1", "randomx0_10"];

% profile_multi_feature(parameters);

% parameters.test_type = "matcutest";

% profile_multi_feature(parameters);

% parameters.multi_feature = ["plain", "noise_1e-1_no_rotation", "noise_1e-2_no_rotation", ...
%         "noise_1e-3_no_rotation", "noise_1e-4_no_rotation", "noise_1e-5_no_rotation", ...
%         "noise_1e-6_no_rotation", "noise_1e-7_no_rotation", "noise_1e-8_no_rotation"];

% % parameters.multi_feature = ["rotation_noisy_1e-7", "rotation_noisy_1e-8"];

% parameters.test_type = "matcutest";

% profile_multi_feature(parameters);

% parameters.test_type = "s2mpj";

% profile_multi_feature(parameters);

% parameters.multi_feature = ["rotation", "rotation_noisy_1e-1",...
%         "rotation_noisy_1e-2", "rotation_noisy_1e-3", "rotation_noisy_1e-4",...
%         "rotation_noisy_1e-5", "rotation_noisy_1e-6", "rotation_noisy_1e-7",...
%         "rotation_noisy_1e-8"];

% parameters.test_type = "matcutest";

% profile_multi_feature(parameters);

% parameters.test_type = "s2mpj";

% profile_multi_feature(parameters);

% parameters.multi_feature = ["rotation_noisy_1e-1", "rotation_noisy_1e-2", "rotation_noisy_1e-3"];

% profile_multi_feature(parameters);

% parameters.test_type = "s2mpj";

% profile_multi_feature(parameters);
