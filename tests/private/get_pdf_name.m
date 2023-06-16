function [pdfname] = get_pdf_name(parameters, i)

prima_list = ["cobyla", "uobyqa", "newuoa", "bobyqa", "lincoa", "mnewuoa_wrapper"];

if parameters.solvers_invoke(i) == "bds"
    pdfname = parameters.Algorithm(i);
elseif parameters.solvers_invoke(i) == "bds_powell"
    powell_factor_stamp = int2str(int32(-log10(parameters.powell_factor(i))));
    pdfname = strcat("GSDS_Powell", "_", powell_factor_stamp);
elseif any(strcmp(prima_list, parameters.solvers_invoke(i)))
    if strcmp(parameters.solvers_invoke(i), "mnewuoa_wrapper")
        pdfname = "mnewuoa";
    else
        pdfname = parameters.solvers_invoke(i);
    end
elseif parameters.solvers_invoke(i) == "matlab_fminsearch"
    pdfname = strcat("fminsearch", "_", "simplex");
elseif parameters.solvers_invoke(i) == "matlab_patternsearch"
    pdfname = strcat("patternsearch", "_", "gps");
elseif parameters.solvers_invoke(i) == "matlab_fminunc"
    pdfname = strcat("fminunc", "_", parameters.fminunc_type);
elseif parameters.solvers_invoke(i) == "CBDS"
    if parameters.polling_outer == "complete"
        polling_outer = "com";
    else
        polling_outer = "opp";
    end
    if parameters.polling_inner == "complete"
        polling_inner = "com";
    else
        polling_inner = "opp";
    end

    pdfname = strcat("CBDS", "_", parameters.nb_tag(i), "_",...
    polling_outer, num2str(parameters.cycling_outer(i)), polling_inner,...
    num2str(parameters.cycling_inner(i)));
elseif parameters.solvers_invoke(i) == "DSPD"
    pdfname = strcat("DSPD", "_", parameters.randomized_strategy(i));
end

end
