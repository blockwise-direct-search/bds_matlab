function [pdfname] = get_pdf_name(parameters, i)
% GET_PDF_NAME gets the part of pdfname of the i-th solver.
prima_list = ["cobyla", "uobyqa", "newuoa", "bobyqa", "lincoa"];

if parameters.solvers_options(i).solver == "bds"
    pdfname = parameters.solvers_options(i).Algorithm;
    pdfname = strcat(pdfname, "_", ...
    parameters.solvers_options(i).sufficient_decrease_factor_level);

elseif parameters.solvers_options(i).solver == "bds_powell"
    powell_factor_stamp = int2str(int32(-log10(parameters.solvers_options(i).powell_factor)));
    pdfname = strcat("CBDS_Powell", "_", powell_factor_stamp);

elseif any(strcmpi(prima_list, parameters.solvers_options(i).solver))
        pdfname = parameters.solvers_options(i).solver;
        if isfield(parameters.solvers_options(i), "version")
            if strcmpi(parameters.solvers_options(i).version, "old")
                pdfname = strcat(parameters.solvers_options(i).solver, "_", "classical");
            end
        end

elseif parameters.solvers_invoke(i) == "wm_newuoa"
    pdfname = parameters.solvers_options(i);

elseif parameters.solvers_invoke(i) == "matlab_fminsearch"
    pdfname = strcat("fminsearch", "_", "simplex");

elseif parameters.solvers_invoke(i) == "matlab_patternsearch"
    pdfname = strcat("patternsearch", "_", "gps");

elseif parameters.solvers_invoke(i) == "matlab_fminunc"
    pdfname = strcat("fminunc", "_", parameters.solvers_options(i).fminunc_type);
    
end

end
