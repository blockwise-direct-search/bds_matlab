function pdfname = get_pdf_name(parameters, i)
% GET_PDF_NAME gets the part of pdfname of the i-th solver.
%

prima_list = ["cobyla", "uobyqa", "newuoa", "bobyqa", "lincoa"];

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
    
    case "dspd"
        pdfname = "dspd";

    case "bds_powell"
        powell_factor_stamp = int2str(int32(-log10(parameters.solvers_options{i}.powell_factor)));
        pdfname = strcat("CBDS_Powell", "_", powell_factor_stamp);

    case "wm_newuoa"
        pdfname = parameters.solvers_options{i}.solver;

    case "nlopt"
        switch parameters.solvers_options{i}.Algorithm
            case "cobyla"
                pdfname = "nlopt_cobyla";
            case "newuoa"
                pdfname = "nlopt_newuoa";
            case "bobyqa"
                pdfname = "nlopt_bobyqa";
        end

    case "fminsearch"
        pdfname = strcat("fminsearch", "_", "simplex");

    case "lam"
        pdfname = "lam";
        if isfield(parameters.solvers_options{i}, "linesearch_type")
            pdfname = strcat(pdfname, "_", ...
                parameters.solvers_options{i}.linesearch_type);
        end

    case "fminunc"
        pdfname = strcat("fminunc", "_", parameters.solvers_options{i}.fminunc_type);
        
    case "patternsearch"
        pdfname = strcat("patternsearch", "_", "gps");
    
    case "bfo_optimize"
        pdfname = "bfo";
end

if any(strcmpi(prima_list, parameters.solvers_options{i}.solver))
    pdfname = parameters.solvers_options{i}.solver;
    if isfield(parameters.solvers_options{i}, "version")
        if strcmpi(parameters.solvers_options{i}.version, "old")
            pdfname = strcat(parameters.solvers_options{i}.solver, "_", "classical");
        end
    end
end

end
