function pdfname = get_pdf_name(parameters, i)
% GET_PDF_NAME gets the part of pdfname of the i-th solver.
%

switch parameters.solvers_options{i}.solver
    case "bds"
        pdfname = upper(parameters.solvers_options{i}.Algorithm);
        if isfield(parameters.solvers_options{i}, "reduction_factor")
            for j = 1:length(parameters.solvers_options{i}.reduction_factor)
                if (-log10(parameters.solvers_options{i}.reduction_factor(j))) < 10
                    pdfname = strcat(pdfname, "_", "0", ...
                        int2str(int32(-log10(parameters.solvers_options{i}.reduction_factor(j)))));
                else
                    if parameters.solvers_options{i}.reduction_factor(j) == eps
                        pdfname = strcat(pdfname, "_", "eps");
                    elseif parameters.solvers_options{i}.reduction_factor(j) == 0
                        pdfname = strcat(pdfname, "_", "00");
                    else
                        pdfname = strcat(pdfname, "_", ...
                            int2str(int32(-log10(parameters.solvers_options{i}.reduction_factor(j)))));
                    end
                end
            end
        end

        if isfield(parameters.solvers_options{i}, "alpha_init_scaling") &&...
                parameters.solvers_options{i}.alpha_init_scaling
            pdfname = strcat(pdfname, "_", "alpha_init_scaling");
        end

        if isfield(parameters.solvers_options{i}, "forcing_function")
            if strcmp(func2str(parameters.solvers_options{i}.forcing_function), func2str(@(x)x.^2))
                pdfname = strcat(pdfname, "_", "quadratic");
            elseif strcmp(func2str(parameters.solvers_options{i}.forcing_function), func2str(@(x)x.^3))
                pdfname = strcat(pdfname, "_", "cubic");
            end
        end

        if isfield(parameters.solvers_options{i}, "forcing_function_type")
            pdfname = strcat(pdfname, "_", parameters.solvers_options{i}.forcing_function_type);
        end

        if isfield(parameters.solvers_options{i}, "shuffling_period")
            if parameters.solvers_options{i}.shuffling_period < 10
                pdfname = strcat(pdfname, "_", "shuffling_period", "_", "0", ...
                    num2str(parameters.solvers_options{i}.shuffling_period));
            else
                pdfname = strcat(pdfname, "_", "shuffling_period", "_", ...
                    num2str(parameters.solvers_options{i}.shuffling_period));
            end
        end

        if isfield(parameters.solvers_options{i}, "cycling_inner")
            pdfname = strcat(pdfname, "_", "cycling_inner", "_", ...
                num2str(parameters.solvers_options{i}.cycling_inner));
        end

        if isfield(parameters.solvers_options{i}, "replacement_delay")
            if parameters.solvers_options{i}.replacement_delay < 10
                pdfname = strcat(pdfname, "_", "replacement_delay", "_", "0", ...
                    num2str(parameters.solvers_options{i}.replacement_delay));
            else
                pdfname = strcat(pdfname, "_", "replacement_delay", "_", ...
                    num2str(parameters.solvers_options{i}.replacement_delay));
            end
        end

        if isfield(parameters.solvers_options{i}, "StepTolerance")
            pdfname = strcat(pdfname, "_", "StepTolerance", "_", ...
                num2str(-log10(parameters.solvers_options{i}.StepTolerance)));
        end

        if isfield(parameters.solvers_options{i}, "MaxFunctionEvaluations_factor")
            pdfname = strcat(pdfname, "_", "MaxFunctionEvaluations_factor", "_", ...
                num2str(parameters.solvers_options{i}.MaxFunctionEvaluations_factor));
        end

        if isfield(parameters.solvers_options{i}, "block_indices_permuted_init")
            pdfname = strcat(pdfname, "_", "block_permuted_init");
        end

        if isfield(parameters.solvers_options{i}, "direction_set_type")
            pdfname = strcat(pdfname, "_", parameters.solvers_options{i}.direction_set_type);
        end

    case "bds_norma"
        pdfname = "bds_norma";

    case "bds_previous"
        pdfname = "bds_previous";

    case "dspd"
        pdfname = "dspd";
        if isfield(parameters.solvers_options{i}, "num_random_vectors")
            if parameters.solvers_options{i}.num_random_vectors < 10
                pdfname = strcat(pdfname, "_", "0", num2str(parameters.solvers_options{i}.num_random_vectors));
            else
                pdfname = strcat(pdfname, "_", num2str(parameters.solvers_options{i}.num_random_vectors));
            end

        end

    case "bds_powell"
        powell_factor_stamp = int2str(int32(-log10(parameters.solvers_options{i}.powell_factor)));
        pdfname = strcat("CBDS_Powell", "_", powell_factor_stamp);

    case "wm_newuoa"
        pdfname = parameters.solvers_options{i}.solver;

    case "nlopt_wrapper"
        switch parameters.solvers_options{i}.Algorithm
            case "cobyla"
                pdfname = "nlopt_cobyla";
            case "newuoa"
                pdfname = "nlopt_newuoa";
            case "bobyqa"
                pdfname = "nlopt_bobyqa";
            case "simplex"
                pdfname = "nlopt_simplex";
        end

    case "fminsearch_wrapper"
        pdfname = strcat("fminsearch", "_", "simplex");

    case "lam"
        pdfname = "lam";
        if isfield(parameters.solvers_options{i}, "linesearch_type")
            pdfname = strcat(pdfname, "_", ...
                parameters.solvers_options{i}.linesearch_type);
        end

    case "fminunc_wrapper"
        pdfname = strcat("fminunc", "_", parameters.solvers_options{i}.fminunc_type);

    case "nomad_wrapper"
        pdfname = "nomad";

    case "patternsearch"
        pdfname = strcat("patternsearch", "_", "gps");

    case "bfo_wrapper"
        pdfname = "bfo";

    case "prima_wrapper"
        pdfname = parameters.solvers_options{i}.Algorithm;
end

end
