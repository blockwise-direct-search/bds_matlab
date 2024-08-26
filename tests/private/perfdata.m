function perfdata(tau, frec, fmin, options_perf)

nprec = length(tau);
prof_output = cell(1, nprec);
format long
for iprec = 1 : nprec
    options_perf.tau = tau(iprec);
    prof_output{iprec} = perfprof(frec, fmin, options_perf);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Appearance of the plots.
%fontsize = 12;
linewidth = 1;
bleu=[0.16, 0.4470, 0.7410];
vert=[0, 0.6, 0];
colors  = {bleu, "k", "b", "r", vert, bleu, "k", "b", "r", vert};
%lines   = {"-", "-.", "--", ":", "-", "-.", "--", ":", "-", "-."};
lines   = {"-", "-", "-", "-", "-", "-", "-", "-", "-", "-"};
hfig = figure("visible", false, "DefaultAxesPosition", [0, 0, 1, 1]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iprec = 1 : nprec
    subplot(1, nprec, iprec);
    if iprec > nprec
        lw = linewidth;
    else
        lw = linewidth * 3;
    end
    for is = 1 : length(options_perf.solvers)
        plot(prof_output{iprec}.profile{is}(1,:), prof_output{iprec}.profile{is}(2,:), ...
            lines{is}, "Color", colors{is},  "Linewidth", lw);
        hold on;
    end
    if iprec == ceil(nprec/2)
        title(strrep(options_perf.feature, '_', '-'), 'FontSize', 14, 'FontWeight', 'bold');
    end
    xlabel(sprintf("%d", iprec));
    axis([0 prof_output{iprec}.cut_ratio 0 1]);
    grid on;
    %pbaspect([1 1 1]);
end
for iprec = 1 : nprec
    ha = get(gcf,"children");
    real_iprec = mod(iprec-1, nprec) + 1;  % The real index of the precision.
    %if iprec > nprec
    %    set(ha(iprec),"position", [0.01+(nprec-real_iprec)/(nprec), 0.53, 0.9/(nprec), 0.4]);
    %else
    %    set(ha(iprec),"position", [0.01+(nprec-real_iprec)/(nprec), 0.1, 0.9/(nprec), 0.4]);
    %end
    if iprec > nprec
        set(ha(iprec),"position", [0.01+(nprec-real_iprec)/(nprec), 0.53, 0.8/(nprec), 0.35]);
    else
        set(ha(iprec),"position", [0.01+(nprec-real_iprec)/(nprec), 0.1, 0.8/(nprec), 0.35]);
    end
end

legend(options_perf.solvers,"Location", "southeast","Orientation","vertical");

% Save the figure as eps.
outdir = fileparts(options_perf.outdir);
figname = strcat("merged", "_", options_perf.pdfname);
epsname = fullfile(outdir, strcat(figname,".eps"));
set(gcf,"position",[0, 0, 4600, 920]);
saveas(hfig, epsname, "epsc2");

% Try converting the eps to pdf.
system(("epstopdf "+epsname));
end

