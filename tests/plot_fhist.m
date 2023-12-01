function plot_fhist(dim, fhist, problem_name, i_run, parameters)
% This file is to draw the function value history of the solvers.
%

color_set = ["red", "blue", "green", "yellow"];
solvers_num = length(parameters.solvers_name);
if iscell(parameters.savepath)
    savepath = parameters.savepath{i_run};
else
    savepath = parameters.savepath;
end
fhist_plot = cell(1, solvers_num);

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

title(char(strcat(num2str(dim), '-', num2str(i_run), "-", problem_name)));
legend(parameters.solvers_name(1), parameters.solvers_name(2), 'Location', 'southwest', 'FontSize', 8);

if dim < 10
    filename = strcat("0", num2str(dim), "_", num2str(i_run), "_", problem_name);
else
    filename = strcat(num2str(dim), "_", num2str(i_run), "_", problem_name);
end

epsname = fullfile(savepath, strcat(filename,'.eps'));
saveas(hfig, epsname, 'epsc2');
% Try converting the eps to pdf.
system(('epstopdf '+epsname));

end

