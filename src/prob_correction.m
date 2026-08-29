function prob = prob_correction(prob,prob_list)

% add model name if necessary
if ~isfield(prob,'modelname')
    prob.modelname = prob_list.name(1:end-4);
end

% remove constraints name to avoid gurobi complaining after removing
% variables or constraints
if isfield(prob,'constrnames')
    prob = rmfield(prob, 'constrnames');
end

if isfield(prob,'varnames')
    prob = rmfield(prob, 'varnames');
end

% there are two bad models in the library
if strcmp(prob.modelname,'mas74')||strcmp(prob.modelname,'mas76')
    % the upper bound of variable x_{151} is 1e12
    % we scale the data for numerical stability
    prob.ub(151) = prob.ub(151)/1e6;
    prob.A(:,151) = prob.A(:,151)*1e6;
end

% some of the data is corrupted, need to fix it as follows
prob.A = sparse(full(prob.A));

% reorder the binary variables to the first
prob = reorder_bvar(prob);

end


function prob = reorder_bvar(prob)
% we reorder the binary variables are arranged as the first k variables
% and the remaining are continuous variables

% number of binary variable
bvar = get_bvar(prob);

% reorder the variables so that the binary variables are the first nb
% variables.
prob.A = [prob.A(:,bvar) prob.A(:,~bvar)];
prob.lb = [prob.lb(bvar); prob.lb(~bvar)];
prob.ub = [prob.ub(bvar); prob.ub(~bvar)];
prob.obj = [prob.obj(bvar); prob.obj(~bvar)];
prob.vtype = [prob.vtype(bvar); prob.vtype(~bvar)];
if isfield(prob,'varnames')
    prob.varnames = [prob.varnames(bvar); prob.varnames(~bvar)];
end

end