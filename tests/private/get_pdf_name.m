function pdfname = get_pdf_name(parameters, i)
% GET_PDF_NAME gets the part of pdfname of the i-th solver.
%

switch parameters.solvers_options{i}.solver
    case "bds"
        pdfname = upper(parameters.solvers_options{i}.Algorithm);
        if isfield(parameters.solvers_options{i}, "sufficient_decrease_factor")
            if parameters.solvers_options{i}.sufficient_decrease_factor == 0
                pdfname = strcat(pdfname, "_", ...
                    num2str(parameters.solvers_options{i}.sufficient_decrease_factor));
            elseif parameters.solvers_options{i}.sufficient_decrease_factor == eps
                pdfname = strcat(pdfname, "_", "eps");
            else
                pdfname = strcat(pdfname, "_", ...
                    int2str(int32(-log10(parameters.solvers_options{i}.sufficient_decrease_factor))));
            end
        end

        if isfield(parameters.solvers_options{i}, "alpha_init_perturbed") &&...
                parameters.solvers_options{i}.alpha_init_perturbed
            pdfname = strcat(pdfname, "_", "perturbed");
        end
    
    case "dspd"
        pdfname = "dspd";

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
