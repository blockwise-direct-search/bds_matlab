function [pdfname] = get_pdf_name(parameters, i)
% GET_PDF_NAME gets the part of pdfname of the i-th solver.
prima_list = ["cobyla", "uobyqa", "newuoa", "bobyqa", "lincoa"];

if parameters.solvers_invoke(i) == "bds"
    pdfname = parameters.Algorithm(i);

elseif parameters.solvers_invoke(i) == "bds_powell"
    powell_factor_stamp = int2str(int32(-log10(parameters.powell_factor(i))));
    pdfname = strcat("CBDS_Powell", "_", powell_factor_stamp);

elseif parameters.solvers_invoke(i) == "bds_cunxin"
    cunxin_factor_tmp = parameters.cunxin_factor{i};
    cunxin_factor_length = length(cunxin_factor_tmp);
    cunxin_factor_stamp = [];

    for j = 1:cunxin_factor_length
        cunxin_factor_stamp = strcat(cunxin_factor_stamp, "_", ... 
            int2str(int32(-log10(cunxin_factor_tmp(j)))));
    end
    pdfname = strcat("CBDS_Cunxin", "_", cunxin_factor_stamp);

elseif any(strcmpi(prima_list, parameters.solvers_invoke(i)))
        pdfname = parameters.solvers_invoke(i);
        if isfield(parameters, "version")
            if strcmpi(parameters.version, "old")
                pdfname = strcat(parameters.solvers_invoke(i), "_", "classical");
            end
        end

elseif parameters.solvers_invoke(i) == "wm_newuoa"
    pdfname = parameters.solvers_invoke(i);

elseif parameters.solvers_invoke(i) == "matlab_fminsearch"
    pdfname = strcat("fminsearch", "_", "simplex");

elseif parameters.solvers_invoke(i) == "matlab_patternsearch"
    pdfname = strcat("patternsearch", "_", "gps");

elseif parameters.solvers_invoke(i) == "matlab_fminunc"
    pdfname = strcat("fminunc", "_", parameters.fminunc_type);
    
end

end
