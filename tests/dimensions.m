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
s.blacklist = [];

% % Problems that take too long to solve.
% % {'FBRAIN3LS'} and {'STRATEC'} take too long for fminunc.
% if ismember("matlab_fminunc", parameters.solvers_name)
%     s.blacklist = [s.blacklist, {'FBRAIN3LS'}, {'STRATEC'}];
% end
% % {"MUONSINELS"} takes nlopt_newuoa so long to run (even making MATLAB crash).
% % {"LRCOVTYPE"}, {'HIMMELBH'} and {'HAIRY'} take nlopt_cobyla so long
% % to run (even making MATLAB crash).
% % {"MUONSINELS"} takes nlopt_bobyqa so long to run (even making MATLAB crash).
% if ismember("nlopt", parameters.solvers_name)
%     s.blacklist = [s.blacklist, {'MUONSINELS'}, {'BENNETT5LS'},...
%         {'HIMMELBH'}, {'HAIRY'}];
% end

if s.mindim >= 6
    s.blacklist = [s.blacklist, { 'ARGTRIGLS', 'BROWNAL', ...
        'COATING', 'DIAMON2DLS', 'DIAMON3DLS', 'DMN15102LS', ...
        'DMN15103LS', 'DMN15332LS', 'DMN15333LS', 'DMN37142LS', ...
        'DMN37143LS', 'ERRINRSM', 'HYDC20LS', 'LRA9A', ...
        'LRCOVTYPE', 'LUKSAN12LS', 'LUKSAN14LS', 'LUKSAN17LS', 'LUKSAN21LS', ...
        'LUKSAN22LS', 'MANCINO', 'PENALTY2', 'PENALTY3', 'VARDIM',
        }];
end

prob_list = secup(s);

end

