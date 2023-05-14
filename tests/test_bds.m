function [] = test_bds(testoptions)

% Get list of problems
s.type = testoptions.problems_type; % Unconstrained: 'u'
s.mindim = testoptions.problems_mindim; % Minimum of dimension
s.maxdim = testoptions.problems_maxdim; % Maximum of dimension
s.blacklist = [];
% Problems that crash.
s.blacklist = [s.blacklist, {}];
% Problems that takes too long to solve.
% {'FBRAIN3LS'} and {'STRATEC'} take too long for fminunc(not for ds and bds).
% {'LRCOVTYPE'} and {'LRIJCNN1'} take long for ds and bds(not for fminunc).
% {'PALMER1C'},{'PALMER2C'},{'PALMER3C'},{'PALMER4C'},{'PALMER5C'},{'PALMER6C'},{'PALMER7C'},{'PALMER8C'}
s.blacklist = [s.blacklist,{'LRCOVTYPE'},{'LRIJCNN1'},{'PARKCH'},{'STRATEC'}...
    ];
list = secup(s);
num = length(list);

for i = 1:num
    p = macup(list(1,i));
    x0 = p.x0;
    dim = length(x0);
    
end


 




end

