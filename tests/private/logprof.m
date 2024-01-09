function logprof(frec, fmin, solvers_name, num_problems, options)
% This file is cited from 
% https://github.com/libprima/prima/blob/main/matlab/tests/private/perfprof.m, which is
% written by Zaikun ZHANG.
%
% This function plots the log-profiles of solvers.
% frec: trajectory of function values; frec(ip, is, ir, k) is the function value of the ip-th
% problem obtained by the is-th solver at the ir-th random run at the k-th iteration.
% fmin: the minimal function values; either fmin(ip) is the minimal function value of the ip-th
% problem, or fmin(ip, ir) is the minimal function value of the ip-th problem for the ir-th run.
% tau: the StepTolerance of convergence.
% solvers: the list of solvers. For logprof, we only compare two solvers.
% options: the options for plotting.

[np, ns, nr, MaxFunctionEvaluations] = size(frec);
% nr is the number of random runs. For logprof, we only compare two solvers and 
% draw the log-profile for the first random run.
nr = 1;
if (ns ~= 2)
    error('logprof only supports two solvers.');
end

% nf_return(ip, is, ir) is the number of function evaluations that the is-th solver uses when it
% returns from solving the ip-th problem at the ir-th random run, and f_return(ip, is, ir) is the
% function value it returns. In testcu.m, the returned function value and constraint violation are
% recorded in fval_history(nf + 1) and cv_history(nf + 1), respectively.
% N.B.: If the problem has no noise, then a reasonable solver (e.g., those in PRIMA) should
% return the best point found along the iterations, in terms of the objective function value or
% a merit function. It is not the case when there is noise.
nf_return = NaN(np, ns, nr);
f_return = NaN(np, ns, nr);
for ip = 1:np
    for is = 1:ns
        for ir = 1:nr
            if all(isnan(frec(ip, is, ir, :)))
                nf_return(ip, is, ir) = MaxFunctionEvaluations;
                f_return(ip, is, ir) = NaN;
            else
                nf_return(ip, is, ir) = find(~isnan(frec(ip, is, ir, :)), 1, 'last') - 1;
                f_return(ip, is, ir) = frec(ip, is, ir, nf_return(ip, is, ir) + 1);
            end
        end
    end
end

% T(ip, is, ir) is the number of function evaluations that the is-th solver needs to solve 
% the ip-th problem (up to StepTolerance tau) at the ir-th random run.
T = NaN(np, ns, nr);
f0 = -Inf(np, nr);
for ip = 1:np
    for ir = 1:nr
        f0(ip,ir) = frec(ip, 1, ir, 1);
    end
end

tau = options.tau;

for ip = 1:np
    for is = 1:ns
        for ir = 1:nr
            if numel(fmin) == length(fmin)  % fmin is a vector indexed by ip only.
                fminp = fmin(ip);
            else
                fminp = fmin(ip, ir);
            end
            fthreshold = tau*f0(ip,ir) + (1-tau)*fminp;
            % We need to ensure fthreshold >= fminp, which may not be true due to rounding
            % errors when fminp = f0(ip,ir).
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            fthreshold = max(fthreshold, fminp);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if options.natural_stop
                % In this case, the number of function evaluations is the amount used by the
                % solver when it stops naturally.
                ftest = f_return(ip, is, ir);
                if (ftest <= fthreshold)
                    T(ip, is, ir) = nf_return(ip, is, ir);
                else
                    T(ip, is, ir) = NaN;
                end
            else
                ftest = min(frec(ip, is, ir, :));
                % Do not change the "if .. else ..." order, as frec(ip, is, ir, 1:M) may be 
                % all NaNs.
                if (ftest <= fthreshold)
                    T(ip, is, ir) = find(frec(ip, is, ir, :) <= fthreshold, 1, 'first');
                else
                    T(ip, is, ir) = MaxFunctionEvaluations;
                end
            end
        end
    end
end

log_ratio = ones(np, nr);
for ip = 1:np
    for ir = 1:nr
        log_ratio(ip, ir) = log2(T(ip, 1, ir)/T(ip, 2, ir));
    end
end
log_ratio = sort(log_ratio);
% Plot the log-profiles.
hfig = figure("visible", "off");
for ir = 1:nr

    bar(log_ratio(:, ir));

    % Set the y-axis range.
    ylim([-20, 20]);

    text(ceil(num_problems/2), 15, char(solvers_name(2)), 'HorizontalAlignment', 'center', 'FontSize', 22);
    text(ceil(num_problems/2), -15, char(solvers_name(1)), 'HorizontalAlignment', 'center', 'FontSize', 22);

    % Add the xlabel.
    xlabel('Problem', 'FontSize', 18);

    % Add the ylabel.
    % Define the original LaTeX string.
    ylabel_prototype = '$log_2(\frac{\mathrm{evals}_{\mathrm{CBDS}}}{\mathrm{evals}_{\mathrm{NEWUOA}}})$';
    % Replace the string using replacement method.
    ylabel_updated = strrep(ylabel_prototype, 'CBDS', char(solvers_name(1)));
    ylabel_updated = strrep(ylabel_updated, 'NEWUOA', char(solvers_name(2)));
    
    % Using the updated LaTeX string to create ylabel.
    ylabel_handle = ylabel(ylabel_updated, 'Interpreter', 'latex');

    % Adjust the position of the ylabel.
    set(ylabel_handle, 'Units', 'normalized', 'Position', [-0.075, 0.5, 0]);

    % Adjust the font size of the ylabel.
    set(ylabel_handle, 'FontSize', 22);

end

% Save the figure as eps.
if int32(-log10(tau)) < 10
    fignamebase = strcat(options.pdfname, '_', 'logperf_', '0', int2str(int32(-log10(tau))));
else
    fignamebase = strcat(options.pdfname, '_', 'logperf_', int2str(int32(-log10(tau))));
end
epsname = fullfile(options.outdir, strcat(fignamebase,'.eps'));
saveas(hfig, epsname, 'epsc2');

% Try converting the eps to pdf.
system(('epstopdf '+epsname));