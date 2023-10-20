function summary_file = perfdata(solvers, frec, fmin, options)


prec = length(options.tau);
nprec = length(prec);
tau = 10.^(-prec);
prof_output = cell(1, 2*nprec);
format long
for iprec = 1 : 2*nprec
    real_iprec = mod(iprec-1, nprec) + 1;  % The real index of the precision.
    prof_options.tau = tau(real_iprec);
    prof_output{iprec} = perfprof(frec, fmin, prof_options);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Appearance of the plots.
%fontsize = 12;
linewidth = 1;
bleu=[0.16, 0.4470, 0.7410];
vert=[0, 0.6, 0];
colors  = {bleu, 'k', 'b', 'r', vert, bleu, 'k', 'b', 'r', vert};
%lines   = {'-', '-.', '--', ':', '-', '-.', '--', ':', '-', '-.'};
lines   = {'-', '-', '-', '-', '-', '-', '-', '-', '-', '-'};
hfig = figure("visible", false, 'DefaultAxesPosition', [0, 0, 1, 1]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iprec = 1 : 2*nprec
    subplot(2, nprec, iprec);
    if iprec > nprec
        lw = linewidth;
    else
        lw = linewidth * 3;
    end
    for is = 1 : ns
       plot(prof_output{iprec}.profile{is}(1,:), prof_output{iprec}.profile{is}(2,:), ...
           lines{is}, 'Color', colors{is},  'Linewidth', lw);
       hold on;
    end
    xlabel(sprintf('%d', iprec));
    axis([0 prof_output{iprec}.cut_ratio 0 1]);
    grid on;
    %pbaspect([1 1 1]);
end
for iprec = 1 : 2*nprec
    ha = get(gcf,'children');
    real_iprec = mod(iprec-1, nprec) + 1;  % The real index of the precision.
    %if iprec > nprec
    %    set(ha(iprec),'position', [0.01+(nprec-real_iprec)/(nprec), 0.53, 0.9/(nprec), 0.4]);
    %else
    %    set(ha(iprec),'position', [0.01+(nprec-real_iprec)/(nprec), 0.1, 0.9/(nprec), 0.4]);
    %end
    if iprec > nprec
        set(ha(iprec),'position', [0.01+(nprec-real_iprec)/(nprec), 0.53, 0.8/(nprec), 0.35]);
    else
        set(ha(iprec),'position', [0.01+(nprec-real_iprec)/(nprec), 0.1, 0.8/(nprec), 0.35]);
    end
end

% The following appears only in the last subplot, but it is sufficient for our use.
ylabel(strrep(test_feature, '_', '\\_'));
legend(solvers,'Location', 'southeast','Orientation','vertical');

% Save the figure as eps.
figname = strcat(stamp, '.', 'summary', '.', feature_and_time);
epsname = fullfile(outdir, strcat(figname,'.eps'));
set(gcf,'position',[0, 0, 4600, 920]);
saveas(hfig, epsname, 'epsc2');

% Try converting the eps to pdf.
try
    system(['epstopdf ',epsname]);
    summary_file = fullfile(outdir, strcat(figname,'.pdf'));
catch
    summary_file = epsname;
end
fprintf('\nSummary for problem type %s with test feature %s:\n\n%s\n\n', options.type, test_feature, summary_file);


% For convenience, save a copy of `problems.txt` and the figures in data_dir. They will be
% replaced in next test with the same `solvers` and `dimrange`.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remove the old files.
delete(fullfile(data_dir, strcat(stamp, '.*.problems.txt')));
delete(fullfile(data_dir, strcat(stamp, '.perf_*.pdf')));
delete(fullfile(data_dir, strcat(stamp, '.perf_*.eps')));
delete(fullfile(data_dir, strcat(stamp, '.summary.*.eps')));
delete(fullfile(data_dir, strcat(stamp, '.summary.*.pdf')));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
copyfile(fprob, data_dir);
epsfiles = dir(fullfile(outdir, '*.eps'));
for k = 1 : length(epsfiles)
    epsname = epsfiles(k).name;
    pdfname = strrep(epsname, '.eps', '.pdf');
    if exist(fullfile(outdir, pdfname), 'file')
        copyfile(fullfile(outdir, pdfname), data_dir);
    else
        copyfile(fullfile(outdir, epsname), data_dir);
    end
end
