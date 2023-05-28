function [pdfname] = get_pdf_name(parameters, i)

prima_list = ["cobyla", "uobyqa", "newuoa", "bobyqa", "lincoa"];

if parameters.solvers_invoke(i) == "blockwise_direct_search"
    pdfname = parameters.solvers_stamp(i);
elseif any(strcmp(prima_list, parameters.solvers_invoke(i)))
    pdfname = parameters.solvers_invoke(i);
elseif parameters.solvers_invoke(i) == "matlab_fminsearch"    
    pdfname = "simplex";
elseif parameters.solvers_invoke(i) == "matlab_fminunc"
    pdfname = strcat("fminunc", "_", parameters.fminunc_type);
elseif parameters.solvers_invoke(i) == "bds_polling"
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
    
    pdfname = strcat("bds_polling", "_", parameters.nb_tag(i), "_",... 
    polling_outer, num2str(parameters.cycling_outer(i)), polling_inner,...
    num2str(parameters.cycling_inner(i)));
end

end

