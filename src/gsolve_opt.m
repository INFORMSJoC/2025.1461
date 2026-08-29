function [x,status_sub] = gsolve_opt(prob,u,d,modelsense,pars)
% optimization version
% status_sub can take the following values
% 'SUCCESS'     we find a new affinely independent solution
% 'INFEASIBLE'  no affinely independent solution exists for this subproblem
% 'TIME_LIMIT'  MILP can't finish within time limit
% 'UNKNOWN'     solves failed for unknown reasons

tol = 1e-2;

% initialize output
% a new affinely indpendent vector if we found it
x = [];

% set up the objective
prob.obj = u;
% adjust the max/min
prob.modelsense = modelsense;
n = length(prob.obj);

% solver setting
gurobi_pars = [];
gurobi_pars.outputflag = 0;
gurobi_pars.TimeLimit = pars.sub_time; % we should never spend too much time in preprocessing
gurobi_pars.PoolSolutions = n;

% solve the optimization problem
output = gurobi(prob,gurobi_pars);

% if the problem is unbounded; then we convert it into a feasibility
% problem
if or(strcmp(output.status,'UNBOUNDED'),strcmp(output.status,'INF_OR_UNBD'))
    prob.obj = zeros(n,1);
    prob.A = [prob.A; u'];
    if strcmp(modelsense,'max')
        prob.rhs = [prob.rhs; d+10];
        prob.sense = [prob.sense; '>'];
    elseif strcmp(modelsense,'min')
        prob.rhs = [prob.rhs; d-10];
        prob.sense = [prob.sense; '<'];
    end
    output = gurobi(prob,gurobi_pars);
end

% if we find at least one solution
if isfield(output,'pool')
    xn = output.pool.xn;
    obj = abs(u'*xn - d);
    [~,idx] = max(obj);
    if obj(idx) > tol
        % remark: of course if we find the optimal solution, then it is
        % always idx
        x = xn(:,idx);
        status_sub = 'SUCCESS';
        return
    end
end

if strcmp(output.status,'TIME_LIMIT')
    status_sub = 'TIME_LIMIT';
    return
end

if strcmp(output.status,'OPTIMAL')
    status_sub = 'INFEASIBLE';
    return
end

status_sub = 'UNKNOWN';

end
