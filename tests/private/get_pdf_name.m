function pdfname = get_pdf_name(parameters, i)
% GET_PDF_NAME gets the part of pdfname of the i-th solver.
%

prima_list = ["cobyla", "uobyqa", "newuoa", "bobyqa", "lincoa"];

switch parameters.solvers_options{i}.solver
    case "bds"
        pdfname = parameters.solvers_options{i}.Algorithm;
        pdfname = strcat(pdfname, "_", ...
            parameters.solvers_options{i}.sufficient_decrease_factor_level);
    case "bds_powell"
        powell_factor_stamp = int2str(int32(-log10(parameters.solvers_options{i}.powell_factor)));
        pdfname = strcat("CBDS_Powell", "_", powell_factor_stamp);
    case "wm_newuoa"
        pdfname = parameters.solvers_options{i}.solver;
    case "nlopt"
        if strcmpi(parameters.solvers_options{i}.Algorithm, "cobyla")
            pdfname = "nlopt_cobyla";
        end
    case "matlab_fminsearch"
        pdfname = strcat("fminsearch", "_", "simplex");
    case "matlab_fminunc"
        pdfname = strcat("fminunc", "_", parameters.solvers_options{i}.fminunc_type);
    case "matlab_patternsearch"
        pdfname = strcat("patternsearch", "_", "gps");
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
