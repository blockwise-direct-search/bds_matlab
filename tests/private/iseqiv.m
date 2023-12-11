function equiv = iseqiv(solvers, p, ir, single_test, prec, options)

pname = p.name;
objective = p.objective;
x0 = p.x0;
n = length(x0);

% Some randomization
% Set seed using pname, n, and ir. We ALTER THE SEED weekly to test the solvers as much as possible.
% N.B.: The weeknum function considers the week containing January 1 to be the first week of the
% year, and increments the number every SUNDAY.
if isfield(options, 'yw')
    yw = options.yw;
elseif isfield(options, 'seed')
    yw = options.seed;
else
    yw = year_week('Asia/Shanghai');
end
fprintf('\nYW = %d\n', yw);
rseed = max(0, min(2^32 - 1,  sum(pname) + n + ir + yw));  % A random seed defined by the current test and yw
orig_rng_state = rng();  % Save the current random number generator settings
rng(rseed);  % Set the random seed for reproducibility
p.x0 = x0 + 0.5*randn(size(x0));
test_options = struct();
if isfield(options, "Algorithm")
    test_options.Algorithm = options.Algorithm;
end
test_options.MaxFunctionEvaluations = max(ceil(20*n*(1+rand)), n+3);  % For reproducibility, do not remove this even if `options` contains `MaxFunctionEvaluations`.
if isfield(options, 'MaxFunctionEvaluations')
    test_options.MaxFunctionEvaluations = options.MaxFunctionEvaluations;
end
test_options.ftarget = objective(x0) - 10*abs(randn)*max(1, objective(x0));
test_options.output_xhist = (rand > 0.5);
test_options.output_block_hist = (rand > 0.5);
if single_test
    % DO NOT INVOKE ANY RANDOMIZATION WITHIN THIS IF. Otherwise, a single test cannot reproduce the
    % corresponding test in a multiple one.
    test_options.output_xhist = true;
    test_options.output_block_hist = true;
end

% Set seed using pname, n, ir, ftarget, and options.Algorithm. We use the seed to make the test reproducible for 
% randomized algorithms (including pbds and rbds).
if isfield(options, 'Algorithm')
    seed = max(0, min(2^32 - 1,  sum(pname) + n + ir + test_options.ftarget + str2double(options.Algorithm) ...
    + norm(p.x0)));
    test_options.seed = seed;
end

if ir == 1
    test_options.expand = exp(abs(randn));
end
if ir == 2
    test_options.shrink = rand;
end
if ir == 3
    test_options.MaxFunctionEvaluations = 1 + ceil(100*n*abs(randn));
end
if ir == 4
    test_options.MaxFunctionEvaluations = 1000*n;
end
if ir == 5
    test_options.MaxFunctionEvaluations = 1;
end
if ir == 6
    test_options.MaxFunctionEvaluations = ceil(n/2);
end
if ir == 7
    test_options.ftarget = inf;
end
if ir == 8
    test_options.alpha_init = 1 + 0.5*(2*rand-1);
    test_options.StepTolerance = test_options.alpha_init;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ir == 9
    test_options.alpha_init = 1 + 0.5*(2*rand-1);
end
if 10 <= ir && ir <= 12
    test_options.StepTolerance = 1e-3*(1 + 0.5*(2*rand-1));
end
if 13 <= ir && ir <= 15
    test_options.reduction_factor = sort(rand(1, 3));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 1 <= ir && ir <= 20
    % The TOUGH tests
    % We must pass the random seed `rseed` to `tough` to ensure reproducibility.
    p = tough(p, rseed);
else
    p.objective  = objective;
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BEGIN: Call the solvers %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% N.B.: In some tests, we may invoke this function with solvers{1} == solvers{2}. So do NOT assume
% that one of the solvers is 'SOLVER' and the other is 'SOLVER_norma'.

% Use function handle to avoid `feval`.
solver1 = str2func(solvers{1});
solver2 = str2func(solvers{2});

test_row_x = (rand > 0.5);
if ~endsWith(solvers{1}, '_norma')
    if test_row_x
        p.x0 = p.x0';
        p.objective = @(x) p.objective(x');
    end
end
if ~endsWith(solvers{2}, '_norma')
    if test_row_x
        p.x0 = p.x0';
        p.objective = @(x) p.objective(x');
    end
end
%tic;
[x1, fx1, exitflag1, output1] = solver1(p.objective, p.x0, test_options);
%T = toc; fprintf('\nRunning time for %s:\t %f\n', solvers{1}, T);
%tic;
[x2, fx2, exitflag2, output2] = solver2(p.objective, p.x0, test_options);
%T = toc; fprintf('\nRunning time for %s:\t %f\n', solvers{2}, T);

% Restore the random number generator state
rng(orig_rng_state);


equiv = iseq(x1(:), fx1, exitflag1, output1, x2(:), fx2, exitflag2, output2, prec);

if ~equiv
    format long;
    fprintf('\nnf: nf1 = %d, nf2 = %d', output1.funcCount, output2.funcCount)
    fprintf('\nx:')
    x1(:)'
    x2(:)'
    (x1(:) == x2(:))'
    fprintf('\nf: fx1 = %.16e, fx2 = %.16e', fx1, fx2)
    fprintf('\nexitflag: exitflag1 = %d, exitflag2 = %d', exitflag1, exitflag2)
    nhist = min(length(output1.fhist), length(output2.fhist));
    fprintf('\nfhist (compare only the last %d evaluations):', nhist);
    output1.fhist
    output2.fhist
    fhist2 = output1.fhist(end-nhist+1: end);
    fhist1 = output2.fhist(end-nhist+1: end);
    fhist1 == fhist2
    if single_test && options.sequential
    %if options.sequential
        fprintf('\nThe solvers produce different results on %s at the %dth run.\n\n', pname, ir);
        cd(options.olddir);
        keyboard
    end
    error('\nThe solvers produce different results on %s at the %dth run.\n', pname, ir);
end

return


function eq = iseq(x, f, exitflag, output, xx, ff, ee, oo, prec)
    eq = true;
    
    if ~isempty(setdiff(fieldnames(output), [fieldnames(oo); 'fhist'; 'xhist'])) ...
            || ~isempty(setdiff(fieldnames(oo), [fieldnames(output); 'fhist'; 'xhist']))
        eq = false;
    end
    
    if (norm(xx-x)/(1+norm(x)) > prec || abs(ff-f)/(1+abs(f)) > prec)
        eq = false;
    end
    
    if isfield(output, 'fhist')
        output.fhist = output.fhist(:);
    else
        output.fhist = [];
    end
    if isfield(oo, 'fhist')
        oo.fhist = oo.fhist(:);
    else
        oo.fhist = [];
    end
    nhist = min(length(output.fhist), length(oo.fhist));
    output.fhist = output.fhist(end - nhist + 1: end);
    oo.fhist = oo.fhist(end - nhist + 1: end);
    
    minfhist = min(length(output.fhist), length(oo.fhist));
    if norm(output.fhist(end-minfhist+1:end) - oo.fhist(end-minfhist+1:end))/(1+norm(output.fhist(end-minfhist+1:end))) > prec
        eq = false;
    end
    
    if (prec == 0 && (exitflag ~= ee|| oo.funcCount ~= output.funcCount))
        eq = false;
    end
    
    return







