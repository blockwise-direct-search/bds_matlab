function [pdfname] = get_pdf_name(parameters, i)

prima_list = ["cobyla", "uobyqa", "newuoa", "bobyqa", "lincoa"];

if parameters.solvers_invoke(i) == "bds"
    pdfname = parameters.Algorithm(i);
elseif parameters.solvers_invoke(i) == "bds_powell"
    powell_factor_stamp = int2str(int32(-log10(parameters.powell_factor(i))));
    pdfname = strcat("GSDS_Powell", "_", powell_factor_stamp);
elseif any(strcmpi(prima_list, parameters.solvers_invoke(i)))
        pdfname = parameters.solvers_invoke(i);
        if isfield(parameters, "version")
            if strcmpi(parameters.version, "old")
                pdfname = strcat(parameters.solvers_invoke(i), "_", "classical");
            end
        end
elseif parameters.solvers_invoke(i) == "matlab_fminsearch"
    pdfname = strcat("fminsearch", "_", "simplex");
elseif parameters.solvers_invoke(i) == "matlab_patternsearch"
    pdfname = strcat("patternsearch", "_", "gps");
elseif parameters.solvers_invoke(i) == "matlab_fminunc"
    pdfname = strcat("fminunc", "_", parameters.fminunc_type);
end

end
