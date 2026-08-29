function bvar = get_bvar(prob)
% get the binary variable indices in prob
bvar = prob.vtype == 'B';

idx = and(prob.lb == 0,prob.ub == 1);
idx = and(idx,prob.vtype == 'I');
bvar(idx) = true; 

end
