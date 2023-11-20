function prob_list = dimensions(type, mindim, maxdim)
fullpath = mfilename("fullpath");
path_examples = fileparts(fullpath);
path_bds = fileparts(path_examples);
path_locate = fullfile(path_bds, "tests", "private");

cd(path_locate)
locate_matcutest();

s.type = type;
s.mindim = mindim;
s.maxdim = maxdim;
prob_list = secup(s);

end

