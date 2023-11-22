function fhist_solvers(problem_name, parameters, test_options)
% This file is to draw the function value history of the solvers.
%

p = macup(problem_name);
x0 = p.x0;
dim = length(x0);

if isfield(parameters, "x0_perturbation_level")
    rr = randn(size(x0));
    rr = rr / norm(rr);
    x0 = x0 + parameters.x0_perturbation_level * max(1, norm(x0)) * rr;
end

color_set = ["red", "blue", "green", "yellow"];
solvers_num = length(parameters.solvers_options);
savepath = parameters.savepath;
fhist = cell(1, solvers_num);
fhist_plot = cell(1, solvers_num);

for i = 1:solvers_num

    solver = str2func(parameters.solvers_options{i}.solver);
    if strcmpi(parameters.solvers_options{i}.solver, "fminunc_wrapper")
        test_options.fd = true;
    else
        test_options.fd = false;
    end
    obj = ScalarFunction(p);
    solver(@(x)obj.fun(x,test_options.is_noisy,1,test_options), x0, parameters.solvers_options{i});
    fhist{i} = obj.valHist;

end

% Get fval.
fval = min(fhist{1});
for i = 1:solvers_num
    if i ~= 1
        fval = min(fval, min(fhist{i}));
    end
end

% Deal fhist.
for i = 1:solvers_num
    fhist_plot{i} = cummin(fhist{i});
end

for i = 1:solvers_num
    fhist_plot{i} = abs(fhist_plot{i} - fval + eps)/max(abs(fhist{i}(1)), eps);
    %fhist_plot{i} = abs(fhist{i} - fval + eps);
end

hfig = figure("visible", false);  % Plot the figure without displaying it.

if parameters.log_x_axis
    for i = 1:solvers_num
        loglog(fhist_plot{i}, color_set(i));
        hold on
    end
else
    for i = 1:solvers_num
        semilogy(fhist_plot{i}, color_set(i));
        hold on
    end
end

title(char(strcat(p.name, '-', num2str(dim))));
legend(parameters.solvers_name(1), parameters.solvers_name(2), 'Location', 'southwest', 'FontSize', 8);

if dim < 10
    filename = strcat("0", num2str(dim), "_", problem_name);
else
    filename = strcat(num2str(dim), "_", problem_name);
end

epsname = fullfile(savepath, strcat(filename,'.eps'));
saveas(hfig, epsname, 'epsc2');

% Try converting the eps to pdf.
system(('epstopdf '+epsname));

end

